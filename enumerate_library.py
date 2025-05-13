#!/usr/bin/env python3
"""
enumerate_library.py  –  N-term/C-term/side-chain enumeration
                         + user‐chosen modes, total‐size preview,
                         sampling, fractional down‐sampling,
                         and accurate progress bar

Generates peptide–organic hybrids ≤150 heavy atoms by attaching motifs
only to the positions you select: N-term, C-term, side-chain, or combinations.
"""

import argparse, csv, gzip, itertools, os, pickle, random, sys
from collections import namedtuple

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from tqdm import tqdm

RDLogger.DisableLog("rdApp.*")


# -----------------------------------------------------------------------------
# Block namedtuple & custom unpickler so pickle.load resolves it
Block = namedtuple("Block", "code smiles description tag")
class BlockUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if name == "Block":
            return Block
        return super().find_class(module, name)


# -----------------------------------------------------------------------------
# Constants & SMARTS
MAX_ATOMS = 150
CHUNK     = 10000
RXN       = AllChem.ReactionFromSmarts("[C:1](=O)O.[N:2]>>[C:1](=O)[N:2]")


# -----------------------------------------------------------------------------
def join_amide(m1, m2):
    prods = RXN.RunReactants((m1, m2))
    if not prods:
        raise ValueError("No reaction product")
    p = prods[0][0]
    Chem.SanitizeMol(p)
    return p

def build_linear(blocks):
    mol = Chem.MolFromSmiles(blocks[0].smiles)
    for blk in blocks[1:]:
        mol = join_amide(mol, Chem.MolFromSmiles(blk.smiles))
    return mol

def load_motifs(pkl_paths):
    motifs = {}
    for fp in pkl_paths:
        with open(fp, "rb") as fh:
            blks = BlockUnpickler(fh).load()
        for b in blks:
            motifs[b.code] = b
    return motifs

HANDLE_MAP = {
    "lys":"Kk","glu":"Ee","asp":"Dd","cys":"Cc",
    "ser":"Ss","thr":"Tt","tyr":"Yy","his":"Hh","trp":"Ww",
}

def sidechain_sites(seq, handles):
    allowed = {h for h in handles if h in HANDLE_MAP}
    return [i for i,res in enumerate(seq)
            if any(res in HANDLE_MAP[h] for h in allowed)]

def embed3d(mol, method):
    """
    Always try to embed at least one conformer;
    then optionally minimize.
    """
    # 1) Embed (ignore ETKDG return code)
    try:
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    except Exception:
        return

    # 2) If embedding failed (no conformer), give up
    if mol.GetNumConformers() == 0:
        return

    # 3) If user asked for no minimization, we're done
    if method == "none":
        return

    # 4) Otherwise try a force-field minimize
    try:
        if method == "mmff":
            props = AllChem.MMFFGetMoleculeProperties(mol)
            ff    = AllChem.MMFFGetMoleculeForceField(mol, props)
        else:
            ff    = AllChem.UFFGetMoleculeForceField(mol)
        if ff is not None:
            ff.Minimize()
    except Exception:
        # best effort; if FF blows up, just leave embedded coords as-is
        return


# -----------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--extra",    nargs="+", required=True,
                  help="pickled blocks (diFunc/monoElec/monoNucl)")
    p.add_argument("--peptide_file",
                  help="one peptide sequence per line")
    p.add_argument("--sample",     action="store_true",
                  help="random sampling mode")
    p.add_argument("--exhaustive", action="store_true",
                  help="full enumeration")
    p.add_argument("--include_d",  action="store_true",
                  help="allow D-amino acids")
    p.add_argument("--min_len", type=int, default=1)
    p.add_argument("--max_len", type=int, default=3)
    p.add_argument("--handles", nargs="*", default=[],
                  help="which sidechains to attach to")
    p.add_argument("--no_motif", action="store_true",
                  help="skip motif attachments")
    p.add_argument("--fraction", type=float, default=1.0,
                  help="down-sample fraction (exhaustive mode)")
    p.add_argument("--estimate",  action="store_true",
                  help="estimate only (exhaustive mode)")
    p.add_argument("--max_enum",  type=int, default=10_000_000,
                  help="cap total ligands (sampling mode)")
    p.add_argument("--minimize", choices=["mmff","uff","none"], default="mmff")
    p.add_argument("--format",   choices=["sdf","pdb","mol2"], default="sdf")
    p.add_argument("--outdir",   default="out")
    p.add_argument("--seed",     type=int, default=42)
    cfg = p.parse_args()

    random.seed(cfg.seed)

    # load AA codes
    from building_blocks import BUILDING_BLOCKS as AA
    alphabet = list(AA.keys()) if cfg.include_d else [c for c in AA if c.isupper()]

    # load motifs
    motifs = load_motifs(cfg.extra)
    di = [m for m in motifs.values() if m.tag=="diFunc"]
    mE = [m for m in motifs.values() if m.tag=="monoElec"]
    mN = [m for m in motifs.values() if m.tag=="monoNucl"]

    # prepare peptide list or generator
    if cfg.peptide_file:
        peptide_list = [list(l.strip()) for l in open(cfg.peptide_file) if l.strip()]
        finite = True
    elif cfg.exhaustive:
        peptide_list = [
            list(prod)
            for length in range(cfg.min_len, cfg.max_len+1)
            for prod in itertools.product(alphabet, repeat=length)
        ]
        finite = True
    else:
        def rand_core():
            while True:
                L = random.randint(cfg.min_len, cfg.max_len)
                yield [random.choice(alphabet) for _ in range(L)]
        peptide_list = None
        cores = rand_core()
        finite = False

    # file‐input summary & mode pick
    if peptide_list is not None:
        Npep = len(peptide_list)
        total_sites = sum(len(sidechain_sites(seq, cfg.handles)) for seq in peptide_list)
        end_motifs  = len(di)+len(mE)+len(mN)
        sc_motifs   = len(mE)+len(mN)

        print(f"\nLoaded {Npep} peptides.")
        print(f"Total handle sites across all peptides: {total_sites}\n")

        print("Select attachment mode:")
        print(" 1) N-terminal only")
        print(" 2) C-terminal only")
        print(" 3) side-chain only")
        print(" 4) N + C terminal only")
        print(" 5) all modes (N, C & side-chain)")
        choice = input("Enter 1–5: ").strip()
        mode = {"1":"N","2":"C","3":"SC","4":"NC","5":"ALL"}.get(choice,"ALL")

        # compute expected total
        if mode in ("N","C"):
            expected = Npep * end_motifs
        elif mode=="NC":
            expected = Npep * 2 * end_motifs
        elif mode=="SC":
            expected = total_sites * sc_motifs
        else:
            expected = Npep * 2 * end_motifs + total_sites * sc_motifs

        print(f"\nExpected total ligands: {expected:,}")
        if input("Proceed? [y/N] ").strip().lower()!="y":
            print("Aborted."); sys.exit(0)
    else:
        mode = "ALL"
        expected = None

    # prepare outputs
    os.makedirs(cfg.outdir, exist_ok=True)
    csv_fh = gzip.open(os.path.join(cfg.outdir,"library.csv.gz"), "wt", newline="")
    csv_w  = csv.writer(csv_fh)
    csv_w.writerow(["sequence","smiles","heavy_atoms"])
    if cfg.format=="sdf":
        sdw = Chem.SDWriter(os.path.join(cfg.outdir,"library.sdf"))
    elif cfg.format=="pdb":
        pdb_dir = os.path.join(cfg.outdir,"pdb"); os.makedirs(pdb_dir, exist_ok=True)

    seen, kept = set(), 0

    # ─── Progress bar setup ───
    # if sampling + peptide file + no custom --max_enum, override max_enum to expected
    default_max = 10_000_000
    if cfg.sample and peptide_list is not None and cfg.max_enum == default_max:
        cfg.max_enum = expected
    bar_total = cfg.max_enum if cfg.sample else expected if peptide_list else None
    pbar = tqdm(total=bar_total, desc="Generating", unit="mol") if bar_total else None

    def write_mol(name, mol):
        nonlocal kept
        mol.SetProp("_Name", name)
        embed3d(mol, cfg.minimize)
        if cfg.format=="sdf":
            sdw.write(mol)
        elif cfg.format=="pdb":
            Chem.MolToPDBFile(mol, os.path.join(pdb_dir,f"{name}.pdb"))
        else:
            Chem.MolToMolFile(mol,
                os.path.join(cfg.outdir,f"{name}.mol2"),
                "mol2", addHs=False, confId=0)
        kept += 1
        if kept % CHUNK == 0:
            csv_fh.flush()
        if pbar:
            pbar.update(1)

    def get_attach_sites(core):
        sites = []
        if mode in ("N","NC","ALL"):
            sites.append(("end",0))
        if mode in ("C","NC","ALL"):
            sites.append(("end",len(core)))
        if mode in ("SC","ALL"):
            for s in sidechain_sites(core, cfg.handles):
                sites.append(("sc",s))
        return sites

    # ─── Sampling mode ───
    if cfg.sample:
        while kept < cfg.max_enum:
            core = random.choice(peptide_list) if peptide_list else next(cores)
            if cfg.no_motif:
                continue
            attach = get_attach_sites(core)
            if not attach:
                continue
            typ,pos = random.choice(attach)
            seq = list(core)
            if typ=="end":
                m = random.choice(di+mE+mN)
                seq.insert(pos, m.code)
            else:
                m = random.choice(mE+mN)
                seq.insert(pos+1, m.code)
            name = "_".join(seq)
            if name in seen:
                continue
            try:
                mol = build_linear([(motifs.get(c) or AA[c]) for c in seq])
            except:
                continue
            if mol.GetNumHeavyAtoms()>MAX_ATOMS:
                continue
            if cfg.fraction<1.0 and random.random()>cfg.fraction:
                continue
            smiles = Chem.MolToSmiles(mol, True)
            csv_w.writerow([name, smiles, mol.GetNumHeavyAtoms()])
            write_mol(name, mol)
            seen.add(name)

    # ─── Exhaustive mode ───
    else:
        def enum_core(core):
            nonlocal kept
            if cfg.no_motif:
                return
            for typ,pos in get_attach_sites(core):
                if typ=="end":
                    for m in (di+mE+mN):
                        seq = list(core); seq.insert(pos, m.code)
                        name="_".join(seq)
                        if name in seen: continue
                        try:
                            mol = build_linear([(motifs.get(c) or AA[c]) for c in seq])
                        except: continue
                        if mol.GetNumHeavyAtoms()>MAX_ATOMS: continue
                        if cfg.fraction<1.0 and random.random()>cfg.fraction: continue
                        csv_w.writerow([name, Chem.MolToSmiles(mol,True), mol.GetNumHeavyAtoms()])
                        write_mol(name, mol); seen.add(name)
                else:
                    for m in (mE+mN):
                        seq = list(core); seq.insert(pos+1, m.code)
                        name="_".join(seq)
                        if name in seen: continue
                        try:
                            mol = build_linear([(motifs.get(c) or AA[c]) for c in seq])
                        except: continue
                        if mol.GetNumHeavyAtoms()>MAX_ATOMS: continue
                        if cfg.fraction<1.0 and random.random()>cfg.fraction: continue
                        csv_w.writerow([name, Chem.MolToSmiles(mol,True), mol.GetNumHeavyAtoms()])
                        write_mol(name, mol); seen.add(name)

        if finite:
            for core in peptide_list:
                if kept>=cfg.max_enum: break
                enum_core(core)
        else:
            while kept < cfg.max_enum:
                enum_core(next(cores))

    # cleanup
    csv_fh.close()
    if cfg.format=="sdf":
        sdw.close()
    if pbar:
        pbar.close()

    print(f"\nFinished – kept {kept:,} molecules   (output → {cfg.outdir})")


if __name__=="__main__":
    main()
