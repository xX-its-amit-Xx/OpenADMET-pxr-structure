#!/usr/bin/env python3
"""
Cross-model GT panel: co-fold the 35 PXR holo crystals (PXR seq + canonical ligand
SMILES) across several OpenProtein models with pooled samples. Scored later vs the
crystal ligand poses -> real per-model accuracy + a test of whether cross-model
features can break the selection wall. Resumable; saves PDBs per model/holo/sample.
"""
import os
import sys
import time
import pandas as pd

OUT = "data/external/posthoc_holo_cofold"
MODELS = ["boltz2", "esmfold2", "af2", "boltz1x"]  # diverse architectures
N_SAMPLES = 5
PXR_SEQ = (
    "GLTEEQRMMIRELMDAQMKTFDTTFSHFKNFRLPGVLSSGCELPESLQAPSREEAAKWSQVRKDLCSLKVSLQLRGEDG"
    "SVWNYKPPADSGGKEIFSLLPHMADMSTYMFKGIISFAKVISYFRDLPIEDQISLLKGAAFELCQLRFNTVFNAETGTW"
    "ECGRLSYCLEDTAGGFQQLLLEPMLKFHYMLKKLQLHEEEYVLMQAISLFSPDRPGVLQHRVVDQLQEQFAITLKSYIE"
    "CNRPQPAHRFLFLKIMAMLTELRSINAQHTQRLLRIQDIHPFATPLMQELFGITGS"
)


def load_holos():
    df = pd.read_csv("data/processed/validation_set/holos_fixed.csv")
    col = "smiles_canonical" if "smiles_canonical" in df else "smiles"
    return {r["holo_id"]: r[col] for _, r in df.iterrows() if isinstance(r[col], str) and r[col]}


def main():
    import openprotein
    from openprotein.molecules.protein import Protein
    from openprotein.molecules.complex import Complex
    from openprotein.molecules.chains import Ligand

    user = os.environ.get("OPENPROTEIN_USERNAME")
    pw = os.environ.get("OPENPROTEIN_PASSWORD")
    if not user or not pw:
        raise SystemExit("Set OPENPROTEIN_USERNAME and OPENPROTEIN_PASSWORD env vars "
                         "(or source a .env) before running.")
    s = openprotein.connect(username=user, password=pw)
    print("connected", flush=True)
    holos = load_holos()
    print(f"{len(holos)} holos with SMILES", flush=True)
    os.makedirs(OUT, exist_ok=True)

    print("building PXR MSA (once)...", flush=True)
    msa = s.align.create_msa(PXR_SEQ)
    msa.wait_until_done()
    print("MSA ready", flush=True)

    for model in MODELS:
        mdir = os.path.join(OUT, model)
        os.makedirs(mdir, exist_ok=True)
        done = {f.split("_s")[0] for f in os.listdir(mdir) if f.endswith(".pdb")}
        todo = {h: sm for h, sm in holos.items() if h not in done}
        print(f"\n=== {model}: {len(done)} done, {len(todo)} todo ===", flush=True)
        fold_fn = getattr(s.fold, model, None)
        if fold_fn is None:
            print(f"  model {model} not available; skip", flush=True)
            continue
        futures = []
        for h, smi in todo.items():
            try:
                prot = Protein.from_expr(PXR_SEQ); prot.set_msa(msa)
                cx = Complex(); cx.set_chain("A", prot); cx.set_chain("X", Ligand(smiles=smi))
                fut = fold_fn.fold(sequences=[cx], diffusion_samples=N_SAMPLES, num_recycles=3)
                futures.append((h, fut))
                print(f"  submitted {h} job={getattr(fut,'job_id','?')}", flush=True)
            except Exception as e:
                print(f"  submit fail {h}: {repr(e)[:120]}", flush=True)
            time.sleep(1)
        # collect
        for h, fut in futures:
            try:
                fut.wait_until_done(timeout=1800)
                res = fut.get()
                # res may be list of structures (one per sample)
                structs = res if isinstance(res, (list, tuple)) else [res]
                for k, st in enumerate(structs):
                    pdb = st.to_string(format="pdb")
                    open(os.path.join(mdir, f"{h}_s{k}.pdb"), "w").write(pdb)
                print(f"  saved {h} ({len(structs)} samples)", flush=True)
            except Exception as e:
                print(f"  collect fail {h}: {repr(e)[:120]}", flush=True)
    print("\nDONE cross-model holo cofold", flush=True)


if __name__ == "__main__":
    main()
