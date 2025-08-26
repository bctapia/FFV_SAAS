"""
ffv.build
Builds the trimer used for all calculations
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Geometry import Point3D

# ---------- tiny helpers ----------

def _pick_label(at: Chem.Atom):
    """Return 'R1'..'R4' if this is a labeled dummy; else None."""
    if at.GetAtomicNum() != 0:
        return None
    # V3000/V2000 R-labels
    for k in ("_MolFileRLabel", "molFileRLabel"):
        if at.HasProp(k):
            try:
                n = at.GetIntProp(k)
            except Exception:
                try: n = int(at.GetProp(k))
                except Exception: n = None
            if n and 1 <= n <= 4:
                return f"R{n}"
    # Pseudo-atom labels
    for k in ("dummyLabel", "molFileAlias"):
        if at.HasProp(k):
            s = at.GetProp(k).strip()
            if s.startswith("R") and s[1:].isdigit():
                n = int(s[1:])
                if 1 <= n <= 4:
                    return f"R{n}"
    return None

def _collect_stubs(mol: Chem.Mol):
    """label -> list of (dummy_idx, neighbor_idx, bondType), sorted by dummy_idx"""
    out = {}
    for at in mol.GetAtoms():
        lab = _pick_label(at)
        if not lab:
            continue
        neighs = [n.GetIdx() for n in at.GetNeighbors()]
        if len(neighs) != 1:
            raise ValueError(f"{lab} must have exactly one neighbor")
        n_idx = neighs[0]
        bt = mol.GetBondBetweenAtoms(at.GetIdx(), n_idx).GetBondType()
        out.setdefault(lab, []).append((at.GetIdx(), n_idx, bt))
    for k in out:
        out[k].sort(key=lambda t: t[0])
    return out

def _ensure_2d(m: Chem.Mol):
    if not m.GetNumConformers():
        AllChem.Compute2DCoords(m)

def _bbox_width(m: Chem.Mol) -> float:
    conf = m.GetConformer()
    xs = [conf.GetAtomPosition(i).x for i in range(m.GetNumAtoms())]
    return (max(xs) - min(xs)) if xs else 0.0

def _translate_inplace(m: Chem.Mol, dx=0.0, dy=0.0, dz=0.0):
    conf = m.GetConformer()
    for i in range(m.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(p.x + dx, p.y + dy, p.z + dz))

def _combine_with_coords(mols):
    """Combine and carry coordinates into a single conformer."""
    combo = mols[0]
    for m in mols[1:]:
        combo = Chem.CombineMols(combo, m)
    counts = [m.GetNumAtoms() for m in mols]
    offs = [0]
    for c in counts[:-1]:
        offs.append(offs[-1] + c)
    conf = Chem.Conformer(sum(counts))
    for off, src in zip(offs, mols):
        sconf = src.GetConformer()
        for i in range(src.GetNumAtoms()):
            conf.SetAtomPosition(off + i, sconf.GetAtomPosition(i))
    out = Chem.Mol(combo)
    out.RemoveAllConformers()
    out.AddConformer(conf)
    return out, offs

# ---------- main: build trimer ----------

def build_trimer_from_repeat_v3000(molblock_v3000: str,
                                   return_mol: bool = False,
                                   add_3d: bool = False,
                                   pad: float = 1.5,
                                   force_single_cross_bonds: bool = False):
    """Auto-stitch a trimer. R1→R2 always; R3→R4 if both exist (ladder)."""
    base = Chem.MolFromMolBlock(molblock_v3000, sanitize=False, removeHs=False)
    if base is None:
        raise ValueError("Could not parse molblock")

    # 2D coords & spacing
    _ensure_2d(base)
    shift = _bbox_width(base) + pad

    # 3 copies in a line
    m1 = Chem.Mol(base)
    m2 = Chem.Mol(base); _translate_inplace(m2, dx=+shift)
    m3 = Chem.Mol(base); _translate_inplace(m3, dx=+2*shift)

    # combine with coordinates
    combo, offs = _combine_with_coords([m1, m2, m3])
    rw = Chem.RWMol(combo)
    off1, off2, off3 = offs[0], offs[1], offs[2]

    # collect stubs on the base once, then shift indices per copy
    stubs = _collect_stubs(base)
    if not stubs.get('R1') or not stubs.get('R2'):
        raise ValueError("Need at least one each of R1 and R2.")
    ladder = bool(stubs.get('R3') and stubs.get('R4'))

    def _shift(lst, off): return [(r+off, n+off, bt) for (r,n,bt) in lst]
    s1 = {k: _shift(v, off1) for k, v in stubs.items()}
    s2 = {k: _shift(v, off2) for k, v in stubs.items()}
    s3 = {k: _shift(v, off3) for k, v in stubs.items()}

    pairs = [('R1','R2')] + ([('R3','R4')] if ladder else [])

    # dedupe edges; safe add
    seen = set()  # set of frozenset({a,b})
    def _safe_add(a, b, bt):
        key = frozenset((a, b))
        if key in seen:
            return
        seen.add(key)
        if rw.GetBondBetweenAtoms(a, b) is None:
            rw.AddBond(a, b, bt)

    def _connect(A, B, la, lb):
        LA, LB = A.get(la, []), B.get(lb, [])
        for (rA, nA, btA), (rB, nB, _btB) in zip(LA, LB):
            bt = Chem.BondType.SINGLE if force_single_cross_bonds else btA
            _safe_add(nA, nB, bt)

    for a, b in pairs:
        _connect(s1, s2, a, b)  # copy1 -> copy2
        _connect(s2, s3, a, b)  # copy2 -> copy3

    # remove dummy R* atoms from all copies
    to_remove = set()
    for S in (s1, s2, s3):
        for lab in ('R1','R2','R3','R4'):
            for (r_idx, _, _) in S.get(lab, []):
                to_remove.add(r_idx)
    for idx in sorted(to_remove, reverse=True):
        rw.RemoveAtom(idx)

    mol = rw.GetMol()
    Chem.SanitizeMol(mol)

    if add_3d:
        try:
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            AllChem.UFFOptimizeMolecule(mol, maxIters=400)
        except Exception:
            pass

    v3k = Chem.MolToMolBlock(mol, forceV3000=True)
    return (mol, v3k) if return_mol else v3k

# ---------- optional: PNG ----------

def save_png(mol: Chem.Mol, path="trimer.png", size=(900, 320), preserve_coords=True):
    if not preserve_coords:
        AllChem.Compute2DCoords(mol)
    Draw.MolToFile(mol, path, size=size)




repeat_mb = """
  -INDIGO-08252523002D

  0  0  0  0  0  0  0  0  0  0  0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 24 24 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 5.2823 -4.60153 0.0 0
M  V30 2 C 7.01296 -4.60078 0.0 0
M  V30 3 C 6.15018 -4.10187 0.0 0
M  V30 4 C 7.01314 -5.60132 0.0 0
M  V30 5 C 5.2833 -5.60499 0.0 0
M  V30 6 C 6.15209 -6.10262 0.0 0
M  V30 7 C 7.88057 -4.10002 0.0 0
M  V30 8 C 4.41819 -4.09997 0.0 0
M  V30 9 R# 8.7466 -4.6003 0.0 0 RGROUPS=(1 2)
M  V30 10 R# 3.55247 -4.59986 0.0 0 RGROUPS=(1 1)
M  V30 11 C 4.41958 -6.10747 0.0 0
M  V30 12 R# 3.55385 -5.60742 0.0 0 RGROUPS=(1 3)
M  V30 13 C 7.88082 -6.09964 0.0 0
M  V30 14 R# 8.74708 -5.59997 0.0 0 RGROUPS=(1 4)
M  V30 15 H 6.15113 -3.10187 0.0 0
M  V30 16 H 6.15561 -7.10262 0.0 0
M  V30 17 H 7.23835 -3.3335 0.0 0
M  V30 18 H 8.52373 -3.33429 0.0 0
M  V30 19 H 5.06237 -3.33509 0.0 0
M  V30 20 H 3.77641 -3.33308 0.0 0
M  V30 21 H 5.064 -6.87214 0.0 0
M  V30 22 H 3.77754 -6.87414 0.0 0
M  V30 23 H 7.2391 -6.86657 0.0 0
M  V30 24 H 8.52398 -6.86537 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 3 1
M  V30 2 2 1 5
M  V30 3 1 5 6
M  V30 4 2 6 4
M  V30 5 1 4 2
M  V30 6 2 2 3
M  V30 7 1 2 7
M  V30 8 1 1 8
M  V30 9 1 7 9
M  V30 10 1 8 10
M  V30 11 1 5 11
M  V30 12 1 11 12
M  V30 13 1 4 13
M  V30 14 1 13 14
M  V30 15 1 3 15
M  V30 16 1 6 16
M  V30 17 1 7 17
M  V30 18 1 7 18
M  V30 19 1 8 19
M  V30 20 1 8 20
M  V30 21 1 11 21
M  V30 22 1 11 22
M  V30 23 1 13 23
M  V30 24 1 13 24
M  V30 END BOND
M  V30 END CTAB
M  END
"""
mol, v3k = build_trimer_from_repeat_v3000(repeat_mb, return_mol=True, add_3d=True)
save_png(mol, "trimer.png")

print(v3k)

#print(trimer_v3000)