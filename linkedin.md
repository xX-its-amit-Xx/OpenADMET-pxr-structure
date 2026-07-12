# LinkedIn post — OpenADMET PXR Structure Challenge wrap-up

> Draft to copy into LinkedIn. Tag the **OpenADMET** and **LatchBio** company pages manually
> for reach (I didn't hardcode @-handles). Links verified live.

---

🧬 **That's a wrap on the OpenADMET PXR Structure Prediction Challenge — I finished rank #2 (LDDT-PLI 0.564, leader 0.573).**

Huge thank you to the **OpenADMET** team for a genuinely hard, genuinely useful challenge: predict the bound pose of 184 ligands (73 PanDDA fragments + 111 drug-like analogs) in the pregnane X receptor (PXR / NR1I2) ligand-binding domain — a notoriously promiscuous, floppy pocket. And a special shoutout to **LatchBio** 💚 for sponsoring the final leg with $500 of compute toward deep AlphaFold3 pose sampling.

Rather than just chase the leaderboard, I treated it as a research problem — and the post-hoc autopsy turned up findings I did *not* expect.

**🎨 The approaches I'm proudest of:**
• **Cross-model z-hybrid selection** — instead of trusting any single model's confidence (pLDDT), I pooled poses across ~15 co-folding models (AF3, Boltz-1/2, OpenFold3, Chai, Protenix, ESMFold2…) and selected on interface-PAE signals. Reframing this from a *generation* problem into a *selection* problem was the single biggest lift.
• **Targeted failure-tail rescue** — my best submission swapped in Protenix-v2 on only the **8 hardest** ligands. Counter-intuitively, *fewer* surgical swaps beat more (8 > 12 > 20).
• A fully open, **per-pose autopsy of all 184 ligands** — with an interactive gallery + 3D viewers.

**🤯 The findings that surprised me most:**
• **The "selection wall."** The pose pool almost always *contained* a near-perfect answer (one holo hid a **0.59 Å** pose) — but no confidence signal, and not even a ground-truth-trained ML selector, could reliably *pick* it. Confidence disagrees with the cross-model consensus in **171 of 184** ligands. It's an information-retrieval failure, not a modeling failure.
• **Input fidelity beats the model zoo.** I co-folded real PXR crystals across three architectures from a *generic* sequence: none placed a single ligand under 2 Å — yet the *same* Boltz with a native, target-specific setup + deep MSA hit 0.59 Å. The lever was never "more models," it was sequence/MSA fidelity and sampling depth.
• **Models cluster by training, not by truth** — architectures agree with each other while all being wrong together. (One caution I only caught in the autopsy: a couple of model exports fold ~20 Å off, which had *inflated* my cross-model disagreement metric ~2× until I re-anchored everything to a real crystal frame. Always check your reference.)

📊 **Full interactive report** (per-pose analysis + the quantified selection wall): https://xx-its-amit-xx.github.io/OpenADMET-pxr-structure/posthoc_analysis.html

🔬 **Per-compound gallery** — all 184, each with a traceable, literature-cited "why we missed" + an interactive 3D cross-model pose scatter: https://xx-its-amit-xx.github.io/OpenADMET-pxr-structure/compound_gallery.html

🎥 **3D pose autopsy** — watch a pose refine toward the crystal (each edit glowing green/red), across three ground-truth holo cases (selection wall vs generation failure): https://xx-its-amit-xx.github.io/OpenADMET-pxr-structure/posthoc_animation.html

🤗 **All 2,246 labeled poses, released CC-BY** so the community can reuse the compute — if your method can pick the good pose already sitting in this pool, that's exactly the open problem: https://huggingface.co/datasets/xX-its-amit-Xx/pxr-structure-pose-pool

💻 **Code:** https://github.com/xX-its-amit-Xx/OpenADMET-pxr-structure

#DrugDiscovery #StructuralBiology #MachineLearning #AlphaFold #Boltz #ComputationalChemistry #OpenScience #PXR #Cheminformatics
