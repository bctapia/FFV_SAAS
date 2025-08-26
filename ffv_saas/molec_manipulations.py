"""
ffv.build
Builds the trimer used for all calculations
"""
from rdkit import Chem
from rdkit.Chem import AllChem
import re

def build_trimer_from_repeat_v3000(molblock_v3000: str, add_3d: bool = False, return_mol: bool = False):
    """
    Build a trimer from a repeat unit drawn with pseudo-atoms R1/R2 (and optionally rails like R1a/R2a, R1b/R2b).
    - The function:
      1) Parses the V3000 repeat unit (keeps pseudo atoms; no sanitize yet)
      2) Identifies R-atom "stubs" and their neighbor + original bond type
      3) Makes 3 copies, stitches: U1.tail -> U2.head, U2.tail -> U3.head (per rail)
      4) Removes all R-atoms
      5) Sanitizes (and optionally embeds + UFF optimizes)
      6) Returns a V3000 Molblock (and optionally the RDKit Mol)

    Heuristics:
      - If both R1 and R2 exist (per rail), R1=head, R2=tail.
      - If only one label exists (e.g., R1 at both ends), the two occurrences are treated as head vs tail by index order.
      - Rails are inferred by optional trailing letters: e.g., R1a/R2a and R1b/R2b form two independent rails.
    """
    # ---- 1) Read (no sanitize yet) ----
    base = Chem.MolFromMolBlock(molblock_v3000, sanitize=False, removeHs=False)
    if base is None:
        raise ValueError("Could not parse V3000 Molblock")

    # ---- helpers ----
    def get_dummy_label(atom):
        if atom.GetAtomicNum() != 0:
            return None
        for key in ("atomLabel", "molFileAlias", "molFileValue", "label", "_MolFileAtomLabel"):
            if atom.HasProp(key):
                v = atom.GetProp(key).strip()
                if v: return v
        # fallback to symbol (often '*')
        return atom.GetSymbol()

    # Map bond types cleanly
    def bond_type_between(m, a1, a2):
        b = m.GetBondBetweenAtoms(a1, a2)
        return b.GetBondType() if b else Chem.BondType.SINGLE

    # Parse label into (base, rail) e.g., "R1a" -> ("R1","a"), "R2" -> ("R2","")
    def split_label(lbl):
        m = re.match(r"^([A-Za-z]+\d*)([A-Za-z]?)$", lbl)  # grab trailing rail letter if present
        if not m: 
            return (lbl, "")
        return (m.group(1), m.group(2))

    # ---- 2) Collect R-stubs: {rail: {label_base: [(r_idx, neigh_idx, bond_type), ...]}} ----
    rails = {}  # rails[""] for single rail; "a","b" for ladders
    for a in base.GetAtoms():
        if a.GetAtomicNum() != 0:
            continue
        lbl = get_dummy_label(a)
        if not lbl:
            continue
        base_label, rail = split_label(lbl)
        # neighbor (should be exactly one)
        neighs = [n.GetIdx() for n in a.GetNeighbors()]
        if len(neighs) != 1:
            # Keep going, but skip malformed stubs
            continue
        stub_idx = a.GetIdx()
        neigh_idx = neighs[0]
        btype = bond_type_between(base, stub_idx, neigh_idx)
        rails.setdefault(rail, {}).setdefault(base_label, []).append((stub_idx, neigh_idx, btype))

    if not rails:
        raise ValueError("No R-labeled pseudo atoms found (R1/R2/etc.)")

    # ---- 3) Make 3 copies and combine ----
    m1 = Chem.Mol(base)
    m2 = Chem.Mol(base)
    m3 = Chem.Mol(base)
    combo = Chem.CombineMols(Chem.CombineMols(m1, m2), m3)
    rw = Chem.RWMol(combo)

    nA = m1.GetNumAtoms()
    nB = m2.GetNumAtoms()
    # Offsets of each unit in the combined graph
    off1, off2, off3 = 0, nA, nA + nB

    # For bond creation we need per-unit stubs with shifted indices
    def shift_stubs(stubs_list, offset):
        # each item: (r_idx, neigh_idx, btype)
        return [(r_idx + offset, neigh_idx + offset, btype) for (r_idx, neigh_idx, btype) in stubs_list]

    # ---- 4) For each rail, determine head/tail and stitch U1->U2 and U2->U3 ----
    # rails[rail] = {"R1":[...], "R2":[...]} or maybe {"R1":[two stubs]} (single label both ends)
    created_bonds = []
    all_dummy_idxs = set()

    for rail, by_label in rails.items():
        labels_present = sorted(by_label.keys())

        # Utility: get stubs for a given label (or synthesize head/tail if only one label present)
        def head_tail_for_unit(unit_by_label):
            """Return (head_stub, tail_stub) where each is a (r_idx, neigh_idx, btype)."""
            # Prefer explicit R1 (head) and R2 (tail)
            if "R1" in unit_by_label and "R2" in unit_by_label:
                head = sorted(unit_by_label["R1"], key=lambda x: x[0])[0]
                tail = sorted(unit_by_label["R2"], key=lambda x: x[0])[-1]
                return head, tail
            # If only R1 exists with two occurrences, use first as head and last as tail
            if "R1" in unit_by_label and len(unit_by_label["R1"]) >= 2:
                ordered = sorted(unit_by_label["R1"], key=lambda x: x[0])
                return ordered[0], ordered[-1]
            # If only R2 exists with two occurrences
            if "R2" in unit_by_label and len(unit_by_label["R2"]) >= 2:
                ordered = sorted(unit_by_label["R2"], key=lambda x: x[0])
                return ordered[0], ordered[-1]
            # Otherwise, pick the lowest-index as head and highest-index as tail from all available
            flat = []
            for L in unit_by_label.values():
                flat.extend(L)
            if len(flat) < 2:
                raise ValueError(f"Rail '{rail}': not enough R-stubs to stitch a chain")
            flat_sorted = sorted(flat, key=lambda x: x[0])
            return flat_sorted[0], flat_sorted[-1]

        # Build unit-local dicts for each copy with shifted indices
        unit1 = {}
        unit2 = {}
        unit3 = {}
        for lbl, lst in by_label.items():
            unit1[lbl] = shift_stubs(lst, off1)
            unit2[lbl] = shift_stubs(lst, off2)
            unit3[lbl] = shift_stubs(lst, off3)

        h1, t1 = head_tail_for_unit(unit1)
        h2, t2 = head_tail_for_unit(unit2)
        h3, t3 = head_tail_for_unit(unit3)

        # Stitch: U1.tail_neighbor -- U2.head_neighbor ; U2.tail_neighbor -- U3.head_neighbor
        # Use the bond type carried by the tail stub (conservative choice)
        _, tail1_nei, tail1_btype = t1
        _, head2_nei, _           = h2
        _, tail2_nei, tail2_btype = t2
        _, head3_nei, _           = h3

        # Create bonds
        rw.AddBond(tail1_nei, head2_nei, tail1_btype); created_bonds.append((tail1_nei, head2_nei))
        rw.AddBond(tail2_nei, head3_nei, tail2_btype); created_bonds.append((tail2_nei, head3_nei))

        # Track all dummy indices to delete later
        for dset in (unit1, unit2, unit3):
            for lst in dset.values():
                for (r_idx, _, _) in lst:
                    all_dummy_idxs.add(r_idx)

    # ---- 5) Delete all R-atoms (descending index order to keep indices valid) ----
    for idx in sorted(all_dummy_idxs, reverse=True):
        rw.RemoveAtom(idx)

    mol_trimer = rw.GetMol()

    # ---- 6) Sanitize (valence, aromaticity, etc.) ----
    Chem.SanitizeMol(mol_trimer)

    # Optional 3D
    if add_3d:
        try:
            AllChem.EmbedMolecule(mol_trimer, AllChem.ETKDGv3())
            AllChem.UFFOptimizeMolecule(mol_trimer, maxIters=400)
        except Exception:
            pass  # embedding can fail for odd cases; keep 2D graph

    # ---- 7) Return V3000 string (and optionally the Mol) ----
    mb = Chem.MolToMolBlock(mol_trimer, forceV3000=True)
    return (mol_trimer, mb) if return_mol else mb


repeat_mb = """
 -INDIGO-08232519232D

  0  0  0  0  0  0  0  0  0  0  0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 14 15 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 7.10634 -7.7625 0.0 0
M  V30 2 C 7.07724 -8.7605 0.0 0
M  V30 3 C 7.92928 -9.28375 0.0 0
M  V30 4 C 8.81045 -8.80901 0.0 0
M  V30 5 C 8.83955 -7.81102 0.0 0
M  V30 6 C 7.9875 -7.28776 0.0 0
M  V30 7 C 9.6625 -9.33226 0.0 0
M  V30 8 C 10.5436 -8.85751 0.0 0
M  V30 9 C 10.5727 -7.85952 0.0 0
M  V30 10 C 9.72072 -7.33627 0.0 0
M  V30 11 R# 6.25429 -7.23925 0.0 0 RGROUPS=(1 1)
M  V30 12 R# 6.19608 -9.23524 0.0 0 RGROUPS=(1 2)
M  V30 13 R# 11.4539 -7.38477 0.0 0 RGROUPS=(1 1)
M  V30 14 R# 11.3956 -9.38075 0.0 0 RGROUPS=(1 2)
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 7 1 4 7
M  V30 8 2 7 8
M  V30 9 1 8 9
M  V30 10 2 9 10
M  V30 11 1 10 5
M  V30 12 1 1 11
M  V30 13 1 2 12
M  V30 14 1 9 13
M  V30 15 1 8 14
M  V30 END BOND
M  V30 END CTAB
M  END
"""
mol_obj, trimer_v3000 = build_trimer_from_repeat_v3000(repeat_mb, return_mol=True)

print(trimer_v3000)