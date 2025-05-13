# building_blocks.py  – 20 L + 20 D α‑amino acids (stereochemistry correct)
from collections import namedtuple
Block = namedtuple("Block", "code smiles description")

aa_L = [
    Block("A", "N[C@@]([H])(C)C(=O)O", "L‑Ala"),
    Block("R", "N[C@@]([H])(CCCNC(=N)N)C(=O)O", "L‑Arg"),
    Block("N", "N[C@@]([H])(CC(=O)N)C(=O)O", "L‑Asn"),
    Block("D", "N[C@@]([H])(CC(=O)O)C(=O)O", "L‑Asp"),
    Block("C", "N[C@@]([H])(CS)C(=O)O", "L‑Cys"),
    Block("E", "N[C@@]([H])(CCC(=O)O)C(=O)O", "L‑Glu"),
    Block("Q", "N[C@@]([H])(CCC(=O)N)C(=O)O", "L‑Gln"),
    Block("G", "NCC(=O)O",                      "Gly"),
    Block("H", "N[C@@]([H])(CC1=CN=CN1)C(=O)O",    "L‑His"),
    Block("I", "N[C@@]([H])(C(C)CC)C(=O)O",        "L‑Ile"),
    Block("L", "N[C@@]([H])(CC(C)C)C(=O)O",        "L‑Leu"),
    Block("K", "N[C@@]([H])(CCCCN)C(=O)O",         "L‑Lys"),
    Block("M", "N[C@@]([H])(CCSC)C(=O)O",          "L‑Met"),
    Block("F", "N[C@@]([H])(Cc1ccccc1)C(=O)O",     "L‑Phe"),
    Block("P", "N1[C@@]([H])(CCC1)C(=O)O",            "L‑Pro"),
    Block("S", "N[C@@]([H])(CO)C(=O)O",            "L‑Ser"),
    Block("T", "N[C@@]([H])(C(O)C)C(=O)O",         "L‑Thr"),
    Block("W", "N[C@@]([H])(CC(=CN2)C1=C2C=CC=C1)C(=O)O", "L‑Trp"),
    Block("Y", "N[C@@]([H])(Cc1ccc(O)cc1)C(=O)O",  "L‑Tyr"),
    Block("V", "N[C@@]([H])(C(C)C)C(=O)O",         "L‑Val"),
]

def invert_alpha(smi: str) -> str:
    if "@@" in smi:  # change L → D
        return smi.replace("@@", "@", 1)
    if "@(" in smi:  # Gly (achiral)
        return smi
    return smi.replace("@", "@@", 1)

aa_D = [Block(b.code.lower(),
              invert_alpha(b.smiles),
              "D‑" + b.description[2:]) for b in aa_L]

BUILDING_BLOCKS = {blk.code: blk for blk in (aa_L + aa_D)}
