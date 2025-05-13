#!/usr/bin/env python3
"""
import_chemdiv_blocks.py
------------------------
Classify a ChemDiv Building‑Blocks SDF into:
    diFunc   – ≥1 electrophile + ≥1 nucleophile (backbone swap capable)
    monoElec – exactly 1 electrophile           (caps N‑terminus / Lys side chain)
    monoNucl – exactly 1 nucleophile            (caps C‑terminus / Glu/Asp side chain)

Output: pickled list of Block(code, smiles, description, tag) used by the enumerator.

Usage example:
python import_chemdiv_blocks.py \
       --in  BBlocks_74721.sdf \
       --out blocks_bb.pkl \
       --max_atoms 150
"""
import argparse, gzip, pickle, collections
from rdkit import Chem

Block = collections.namedtuple("Block", "code smiles description tag")

SMARTS = {
    "carboxyl":   "[CX3](=O)[OX1H0-,OX2H1]",
    "amine":      "[NX3;H2,H1;!$(NC=O)]",
    "isocyanate": "N=C=O",
    "chloroform": "COC(=O)Cl",
    "sulfonylCl": "[SX4](=O)(=O)Cl",
}

def detect_tag(mol):
    hit = {k: len(mol.GetSubstructMatches(Chem.MolFromSmarts(s))) for k, s in SMARTS.items()}
    elec = hit["carboxyl"] + hit["isocyanate"] + hit["chloroform"] + hit["sulfonylCl"]
    nucl = hit["amine"]
    if elec >= 1 and nucl >= 1:
        return "diFunc"
    if elec == 1 and nucl == 0:
        return "monoElec"
    if nucl == 1 and elec == 0:
        return "monoNucl"
    return None

def mol_supplier(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rb") as fh:          # binary mode fix
        for mol in Chem.ForwardSDMolSupplier(fh):
            if mol:
                yield mol

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in",  dest="sdf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max_atoms", type=int, default=150)
    args = ap.parse_args()

    stats = collections.Counter()
    blocks = []
    for idx, mol in enumerate(mol_supplier(args.sdf), start=1):
        if mol.GetNumHeavyAtoms() > args.max_atoms:
            stats["too_big"] += 1
            continue
        tag = detect_tag(mol)
        if tag is None:
            stats["no_handle"] += 1
            continue
        code = f"X{idx:06d}"
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else code
        blocks.append(Block(code, Chem.MolToSmiles(mol, True), name + ":" + tag, tag))
        stats[tag] += 1

    print("=== import summary ===")
    for k, v in stats.items():
        print(f"{k:12s}: {v}")
    print("kept total:", len(blocks))

    with open(args.out, "wb") as fh:
        pickle.dump(blocks, fh)

if __name__ == "__main__":
    main()
