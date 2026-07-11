# Why we missed — PXR pose-failure research provenance

Deep-research workflow (6 angles, 25 sources, 106 claims, 3-vote adversarial verification -> 24 confirmed, 1 refuted). Backs the per-site hover explanations in docs/posthoc_animation.html.

## Verified claims (2/3+ votes)

1. The PXR ligand-binding pocket is large (~1,200 A^3) and expandable to >1,500 A^3, molding itself to diverse ligand shapes — a physical basis for near-isoenergetic alternate poses that defeat pose prediction.  
   _source: https://www.cell.com/structure/fulltext/S0969-2126(23)00297-6_
2. PXR selectivity/binding arises jointly from binding-pocket malleability and ligand conformational flexibility, meaning a single ligand can adopt genuinely different low-energy binding modes in the same pocket.  
   _source: https://www.cell.com/structure/fulltext/S0969-2126(23)00297-6_
3. His407 acts as a polar anchor in the PXR pocket, accepting a hydrogen bond from a ligand hydroxyl/polar group — directly supporting the site-2 His407 anchoring role in the research question.  
   _source: https://www.cell.com/structure/fulltext/S0969-2126(23)00297-6_
4. The PXR ligand-binding cavity is buried with a volume greater than 1,000 A^3, notably larger than the pockets of other nuclear receptors, directly supporting the 'large ~1150 A^3 promiscuous flexible pocket' premise.  
   _source: https://www.nature.com/articles/s41598-018-34373-z_
5. The PXR pocket interior is essentially uncharged and hydrophobic (20 of the cavity-lining residues hydrophobic, only 4 polar and 4 charged, with an E321-R410 salt bridge neutralizing charge), explaining weak directional anchoring and hydrophobic sub-pocket degeneracy that misleads co-folding models.  
   _source: https://www.nature.com/articles/s41598-018-34373-z_
6. SR12813, a canonical PXR agonist, anchors via a polar cluster comprising S247, Q285, and H407 (plus W299/R410) forming hydrogen bonds to the ligand's hydroxyl/phosphate groups, confirming these exact residues as the PXR polar anchor cluster the models displaced.  
   _source: https://www.nature.com/articles/s41598-018-34373-z_
7. A single PXR ligand can adopt multiple near-isoenergetic binding geometries: SR12813 was observed in five different orientations within the cavity across experimental crystal structures, evidencing the pose degeneracy that defeats single-answer pose prediction.  
   _source: https://www.nature.com/articles/s41598-018-34373-z_
8. The PXR ligand-binding pocket is unusually large (~1,200 Å³) and can expand to >1,500 Å³ to accommodate large or diverse ligands, which underlies its promiscuity and makes a single correct pose hard to pin down.  
   _source: https://www.sciencedirect.com/science/article/pii/S0969212623002976_
9. PXR's large pocket (1200-1600 Å³) permits binding of a vast structural array of molecules, i.e. it is a promiscuous, low-specificity pocket rather than a lock-and-key site.  
   _source: https://www.sciencedirect.com/science/article/pii/S0969212623002976_
10. Pocket malleability combined with ligand flexibility lets PXR mold its pocket to diverse ligand shapes and support entirely new binding modes, meaning multiple near-equivalent poses are geometrically feasible.  
   _source: https://www.sciencedirect.com/science/article/pii/S0969212623002976_
11. Binding of hyperforin induces conformational changes in the PXR ligand-binding pocket relative to prior human PXR structures and expands the pocket volume by 250 Å^3, demonstrating the pocket is not rigid but reshapes around different ligands.  
   _source: https://pubs.acs.org/doi/10.1021/bi0268753_
12. PXR uses structural flexibility of its ligand-binding domain as its core mechanism for recognizing a broad, structurally diverse chemical space of ligands.  
   _source: https://pubs.acs.org/doi/10.1021/bi0268753_
13. AF3's dominant docking failure mode is ligand-pocket orientation error (misplacing the ligand within a correctly identified pocket), accounting for 58.6% of errors on unseen ligands — far more than ligand-structure (29.9%) or pocket-structure (17.2%) errors.  
   _source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12424161/_
14. AF3 is biased toward placing ligands so that headgroups interact with polar/charged residues and misplaces hydrophobic moieties, misaligning the hydrophobic tail relative to the native binding mode.  
   _source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12424161/_
15. The PXR ligand-binding cavity is a large ~1150 Å³ pocket in the apo state that expands to accommodate larger ligands (up to ~1544 Å³ with hyperforin), giving it 'directed promiscuity' — the ability to bind structurally diverse ligands.  
   _source: https://pmc.ncbi.nlm.nih.gov/articles/PMC2789303/_
16. Hot-spot/probe hydrogen-bonding analysis identifies Gln285, Ser247, His407, His327, and Cys284 as the polar residues mediating the largest fractions of hydrogen bonds, quantitatively establishing them as the pocket's directional anchor cluster.  
   _source: https://pmc.ncbi.nlm.nih.gov/articles/PMC2789303/_
17. Computational mapping of the PXR pocket revealed five distinct binding subsites (Right, Left, Up, Down, Center) into which ligands extend depending on size and shape, explaining broad promiscuity with near-isoenergetic multi-subsite binding modes.  
   _source: https://pmc.ncbi.nlm.nih.gov/articles/PMC2789303/_
18. AI co-folding/docking methods degrade sharply when the target pocket is dissimilar from training pockets — Chai-1's median RMSD rises from ~2.37 Å to ~5.69 Å as pocket similarity drops — evidence of overfitting to training pocket conformations rather than generalizing to novel binding sites.  
   _source: https://arxiv.org/html/2505.01700v2_
19. The PXR ligand-binding pocket is large and conformable, expanding from ~1280 Å³ (SR12813 complex) to more than 1600 Å³, allowing it to accommodate ligands ranging from 250 to >800 Da — explaining the near-isoenergetic, weakly-anchored poses that make co-folding placement ambiguous.  
   _source: https://academic.oup.com/mend/article/19/12/2891/2737785_
20. Only six side chains are consistently involved in ligand binding across all PXR LBD structures: three polar (Ser-247, Gln-285, His-407) and three hydrophobic (Met-243, Trp-299, Phe-420) — pinpointing the exact polar anchor cluster (Ser247/Gln285/His407) the model displaced and the hydrophobic residues (Phe420) lining the sub-pocket.  
   _source: https://academic.oup.com/mend/article/19/12/2891/2737785_
21. PXR exhibits binding promiscuity unlike any other nuclear receptor, binding agonists that vary widely in shape and chemical features and in mass from 250 to >800 Da — the basis for why weak directional anchoring produces multiple competing poses.  
   _source: https://academic.oup.com/mend/article/19/12/2891/2737785_
22. The amide-N-H···S=C hydrogen bond (thiocarbonyl sulfur as acceptor) is conventionally considered weak because of the small electronegativity of sulfur — establishing why C=S is treated as a weak/atypical H-bond acceptor relative to carbonyl oxygen.  
   _source: https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.7b01810_
23. The electrostatic potential around sulfur has both positive (sigma-hole) and negative (lone-pair) regions, so sulfur accepts hydrogen bonds only in specific lone-pair directions rather than isotropically.  
   _source: https://pmc.ncbi.nlm.nih.gov/articles/PMC5816952/_
24. Thiocarbonyl sulfur (C=S) is a weaker hydrogen bond acceptor than carbonyl oxygen (C=O): thioacetone shows an interaction energy of -8.4 kcal/mol with HF versus -10.0 kcal/mol for acetone (~20% weaker).  
   _source: https://pmc.ncbi.nlm.nih.gov/articles/PMC7325729/_

## Refuted (NOT asserted)

- The PXR-agonist SR12813 was co-crystallized bound within the human PXR (NR1I2) ligand-binding domain in complex with RXRalpha, providing the primary structural evidence cited for SR12813's binding mode in the PXR LBD.

## Sources

- https://www.rcsb.org/structure/4J5X
- https://www.cell.com/structure/fulltext/S0969-2126(23)00297-6
- https://www.nature.com/articles/s41598-018-34373-z
- https://www.sciencedirect.com/science/article/pii/S0969212623002976
- https://pubs.acs.org/doi/10.1021/bi0268753
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12424161/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC2789303/
- https://arxiv.org/html/2505.01700v2
- https://academic.oup.com/mend/article/19/12/2891/2737785
- https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.7b01810
- https://pmc.ncbi.nlm.nih.gov/articles/PMC5816952/
- https://par.nsf.gov/servlets/purl/10318435
- https://pmc.ncbi.nlm.nih.gov/articles/PMC7325729/
- https://pubmed.ncbi.nlm.nih.gov/28876948/
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8657569/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC8864553/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC9563780/
- https://academic.oup.com/mend/article/19/5/1125/2737859
- https://www.biorxiv.org/content/10.1101/792671v3.full
- https://febs.onlinelibrary.wiley.com/doi/full/10.1046/j.1432-1033.2002.03207.x
- https://pubs.acs.org/doi/10.1021/ci049722q
- https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1000594
- https://pmc.ncbi.nlm.nih.gov/articles/PMC8390552/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC2781111/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC10872772/