from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd


def mol2xyz(mol, comment=None):
    c = mol.GetConformers()[0]
    coords = c.GetPositions()
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]

    xyz = "{}\n{}\n".format(len(atoms), comment)
    for a, c in zip(atoms, coords):
        xyz += "{0}     {1:14.9f}    {2:14.9f}    {3:14.9f}\n".format(a, *c)

    return xyz


def xyz2mol(xyz, smiles):
    lines = xyz.splitlines()
    N_atoms = int(lines[0])
    comments = lines[1]

    if N_atoms != len(lines[2:]):
        raise ValueError("Number of atoms does not match")

    mol = Chem.MolFromSmiles(smiles)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    mol = Chem.AddHs(mol, addCoords=True)
    try:
        conf = mol.GetConformers()[0]
    except:
        AllChem.EmbedMultipleConfs(
            mol,
            numConfs=1,
            pruneRmsThresh=0.5,
            randomSeed=1,
            useExpTorsionAnglePrefs=True,
            useBasicKnowledge=True,
        )
        try:
            conf = mol.GetConformers()[0]
        except:
            return None, None

    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    for i, coord in enumerate(lines[2:]):
        coord = coord.split()

        if atoms[i] != coord[0]:
            raise ValueError("Atom does not match")

        conf.SetAtomPosition(i, np.array(coord[1:]).astype("float"))

    mol.SetProp("comments", comments)
    return mol, comments


def xyz2com(xyz, head, footer, comfile, charge=0, mult=1, title="Title"):
    coords = [x for x in xyz.splitlines()]
    new_coords = []

    for coord in coords:
        symbol, x, y, z = coord.split()
        new_coords.append("{} {:.6f} {:.6f} {:.6f}\n".format(symbol, float(x), float(y), float(z)))

    with open(comfile, "w") as com:
        com.write(head)
        com.write("\n")
        com.write(title + "\n")
        com.write("\n")
        com.write("{} {}\n".format(charge, mult))
        com.writelines(new_coords)
        com.write("\n")
        com.write(footer)
        com.write("\n\n\n")


def write_mol_to_sdf(mol, path, confIds=[0], confEns=None):
    if isinstance(confIds, int):
        confIds = [confIds]
    if isinstance(confEns, int):
        confEns = [confEns]
    writer = Chem.SDWriter(path)
    if confEns:
        for confId, confEn in zip(confIds, confEns):
            mol.SetProp("ConfEnergies", str(confEn))
            writer.write(mol, confId=confId)
    else:
        for confId in confIds:
            writer.write(mol, confId=confId)
    writer.close()


def write_mols_to_sdf(mols, path):
    writer = Chem.SDWriter(path)
    for mol in mols:
        writer.write(mol, confId=0)
    writer.close()


def load_sdf(path, removeHs=False, sanitize=False):
    return Chem.SDMolSupplier(path, removeHs=removeHs, sanitize=sanitize)

def clean_xyz_str(xyz_str):
    coords = [x for x in xyz_str.splitlines()]
    new_coords = []

    for coord in coords:
        symbol, x, y, z = coord.split()
        new_coords.append("{} {:.6f} {:.6f} {:.6f}\n".format(symbol, float(x), float(y), float(z)))

    return "".join(new_coords)