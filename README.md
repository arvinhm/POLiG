# Peptide‑Organic Library Generator

A turnkey toolkit that builds large‑scale libraries of peptide–organic hybrids (or pure peptides) ready for docking/MD.

---

## Folder structure

```
project/
├── building_blocks.py           # 20 L + 20 D α‑AAs (stereochemistry correct)
├── import_chemdiv_blocks.py     # parses & classifies ChemDiv Building‑Blocks SDF
├── enumerate_library.py         # generates library, SMILES & 3‑D
├── blocks_bb.pkl                # ← created by importer (pickled motifs)
└── lib_*/                       # output directories (one per run)
    ├── library.csv.gz           # sequence,SMILES,heavy_atoms
    └── library.sdf / .pdb /.mol2
```

---

## 1  Set‑up

```bash
conda create -n pep_enum python=3.11 rdkit -c conda-forge -y
conda activate pep_enum
```

---

## 2  Import ChemDiv motifs

```bash
python import_chemdiv_blocks.py \
       --in  ./BBlocks_74721.sdf \
       --out blocks_bb.pkl \
       --max_atoms 75
```

* Prints a handle summary.
* Creates `blocks_bb.pkl` used by the enumerator.

---

## 3  Generate a library

### Truly random peptides

Total number of 42247 organic compund we have ==> 42247 * 20^N will be the number of combinations
if you want to add compund in backbone use --with_backbone_swap flag

```bash
python enumerate_library.py \
  --extra blocks_bb.pkl \
  --min_len 1 --max_len 3 \
  --handles lys glu asp cys ser thr tyr his trp \
  --sample \
  --max_enum 100 \
  --format pdb \
  --outdir lib_sample100

```

### If you want to take 10 % of the full enumeration

```bash
python enumerate_library.py \
  --extra blocks_bb.pkl \
  --min_len 1 --max_len 3 \
  --handles lys glu asp cys ser thr tyr his trp \
  --with_backbone_swap \
  --exhaustive \
  --fraction 0.10 \
  --format pdb \
  --outdir lib_10pct \
  --estimate

```


### C) User‑supplied cores

`peptides.txt` (one sequence per line):

```
MCESFE
ALAILE
dKGV
```
```bash
python enumerate_library.py \
    --peptide_file peptides.txt \
  --extra blocks_bb.pkl \
  --min_len 1 --max_len 3 \
  --handles lys glu asp cys ser thr tyr his trp \
  --format pdb \
  --outdir lib_sample100

```


---
## Peptide lib names

```bash
N_t_S_X004920.pdb
│ │ │ └─ X004920  →  the organic motif code (pulled from blocks_bb.pkl)
│ │ └─── “S”      →  L‑Ser side‑chain that carries the motif, or (if no side‑chain)
│ └──── “t”       →  **D‑Thr**  (lower‑case = D‑isomer, upper‑case = L‑isomer)
└───── “N”        →  L‑Asn      (first residue in the linear peptide core)
```

## Key flags

| Flag                   | Meaning                                             |        |                               |
| ---------------------- | --------------------------------------------------- | ------ | ----------------------------- |
| `--peptide_file FILE`  | Use exactly those cores.                            |        |                               |
| `--exhaustive`         | Full Cartesian enumeration (no random).             |        |                               |
| `--min_len/--max_len`  | Core length range when generating.                  |        |                               |
| `--handles lys glu …`  | Residues allowed for side‑chain attachment.         |        |                               |
| `--with_backbone_swap` | Allow motif to replace one backbone bond.           |        |                               |
| `--no_motif`           | Build peptide‑only library.                         |        |                               |
| `--max_enum N`         | Stop after N final molecules.                       |        |                               |
| `--sample`             | Shuffle cores before enumeration (exhaustive mode). |        |                               |
| \`--format sdf         | pdb                                                 | mol2\` | 3‑D file type.                |
| \`--minimize mmff      | uff                                                 | none\` | Geometry optimisation method. |

---

## Output

* **CSV** – `sequence,smiles,heavy_atoms` (gzip‑compressed).
* **3‑D file** – minimised structures (one record per molecule).
* All molecules are ≤ 150 heavy atoms and RDKit‑sanitised.

---

## Troubleshooting

* **Importer fails** → check SDF path or RDKit install.
* **Enumeration slow** → lower `--max_enum`, use more CPUs, or restrict handles.
* **3‑D write exceeds disk** → lower `--max_out`.

---

## License

