#!/usr/bin/env python
# coding: utf-8

import os
import re
import tarfile
import numpy as np
import rdkit

from rdmc.mol import RDKitMol

from .utils import make_xyz_str

periodictable = [
    "",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Uub",
    "Uut",
    "Uuq",
    "Uup",
    "Uuh",
    "Uus",
    "Uuo",
]


def check_job_status(member, tar):
    f = tar.extractfile(member)
    lines = f.readlines()
    for line in reversed(lines):

        if b"Normal termination" in line:
            return True
        else:
            return False


# In[7]:


def get_cpu(member, tar):
    f = tar.extractfile(member)
    lines = f.readlines()
    for line in lines:
        if line.find(b"Job cpu time") > -1:
            days = int(line.split()[3])
            hours = int(line.split()[5])
            mins = int(line.split()[7])
            secs = float(line.split()[9])
            CPU = tuple([days, hours, mins, secs])
            return CPU


# In[8]:


def get_wall(member, tar):
    f = tar.extractfile(member)
    lines = f.readlines()
    for line in lines:
        if line.find(b"Elapsed time") > -1:
            days = int(line.split()[2])
            hours = int(line.split()[4])
            mins = int(line.split()[6])
            secs = float(line.split()[8])
            CPU = tuple([days, hours, mins, secs])
            return CPU


# In[9]:


def make_input_file_from_xyz(symbols, coords):
    xyz_str = ""
    for s, c in zip(symbols, coords):
        xyz_str = xyz_str + f"{s}  {c[0]: .10f}  {c[1]: .10f}  {c[2]: .10f}\n"
    return xyz_str


# In[10]:


def load_geometry(member, tar, periodictable=periodictable, initial=False):
    """
    Return the optimum geometry of the molecular configuration from the
    Gaussian log file. If multiple such geometries are identified, only the
    last is returned.
    """
    step = -1
    idx, number, coord, symbol = [], [], [], []
    f = tar.extractfile(member)
    line = f.readline()
    while line != b"":
        # Automatically determine the number of atoms
        if b"Input orientation:" in line:
            step += 1
            number, coord = [], []
            for i in range(5):
                line = f.readline()
            while (
                b"---------------------------------------------------------------------"
                not in line
            ):
                data = line.split()
                idx.append(int(data[0]))
                number.append(int(data[1]))
                coord.append([float(data[3]), float(data[4]), float(data[5])])
                line = f.readline()
        line = f.readline()

        if coord and initial:
            break

    number = np.array(number)
    symbol = [periodictable[x] for x in number]

    xyz_dict = dict()
    for x in zip(idx, symbol, coord):
        xyz_dict[x[0]] = (x[1], tuple(x[2]))

    xyz_str = make_xyz_str(symbol, coord)
    return xyz_str, xyz_dict, step


def load_geometry_std(member, tar, periodictable=periodictable, initial=False):
    """
    Return the optimum geometry IN STANDARD ORIENTATION of the molecular configuration from the
    Gaussian log file. If multiple such geometries are identified, only the
    last is returned.
    """
    step = -1
    idx, number, coord, symbol = [], [], [], []
    f = tar.extractfile(member)
    line = f.readline()
    while line != b"":
        # Automatically determine the number of atoms
        if b"Standard orientation:" in line:
            step += 1
            number, coord = [], []
            for i in range(5):
                line = f.readline()
            while (
                b"---------------------------------------------------------------------"
                not in line
            ):
                data = line.split()
                idx.append(int(data[0]))
                number.append(int(data[1]))
                coord.append([float(data[3]), float(data[4]), float(data[5])])
                line = f.readline()
        line = f.readline()

        if coord and initial:
            break

    number = np.array(number)
    symbol = [periodictable[x] for x in number]

    xyz_dict = dict()
    for x in zip(idx, symbol, coord):
        xyz_dict[x[0]] = (x[1], tuple(x[2]))

    xyz_str = make_input_file_from_xyz(symbol, coord)
    return xyz_str, xyz_dict, step


# In[11]:


def load_freq(member, tar):
    """
    Return the frequencies from a transition state frequency
    calculation in cm^-1.
    """
    frequencies = []
    f = tar.extractfile(member)
    line = f.readline()
    while line != b"":
        # Read vibrational frequencies
        if b"Frequencies --" in line:
            frequencies.extend(line.split()[2:])
        line = f.readline()

    frequencies = [float(freq) for freq in frequencies]
    frequencies.sort()

    return frequencies


def load_first_mode(member, tar, periodictable=periodictable):
    f = tar.extractfile(member)
    line = f.readline()

    idx, number, coord, symbol = [], [], [], []

    while line != b"":
        if b"Frequencies --" in line:
            for i in range(5):
                line = f.readline()
            while len(line.split()) > 3:
                data = line.split()
                idx.append(int(data[0]))
                number.append(int(data[1]))
                coord.append([float(data[2]), float(data[3]), float(data[4])])
                line = f.readline()
            else:
                break
        line = f.readline()

        number = np.array(number)
        symbol = [periodictable[x] for x in number]

        result = dict()
        for x in zip(idx, symbol, coord):
            result[x[0]] = (x[1], tuple(x[2]))
        return result


# In[12]:


def check_neg_freq(frequencies):
    neg_idx = np.where(np.array(frequencies) < 0)[0]
    if len(neg_idx) >= 1:
        raise ValueError("Imaginary frequency found")
    else:
        return frequencies


# In[14]:


def check_freq(member, tar):

    freq = load_freq(member, tar)

    try:
        check_neg_freq(freq)
    except:
        return False

    return True


# In[15]:


def load_zpe_and_scf(member, tar):
    f = tar.extractfile(member)
    s = f.read()
    s = s.replace(b"\n", b"").replace(b" ", b"")
    zpe = float(re.findall(b"ZeroPoint=(-*\d+.\d+)", s)[0])
    scf = float(re.findall(b"HF=(-*\d+.\d+)", s)[0])
    return zpe, scf


# In[16]:


def load_e0_zpe(member, tar):
    f = tar.extractfile(member)
    line = f.readline()
    while line != b"":
        if b"Sum of electronic and zero-point Energies=" in line:
            e0_zpe = float(line.split()[-1])
            break
        line = f.readline()
    return float(e0_zpe)


# In[17]:


def load_gibbs(member, tar):
    f = tar.extractfile(member)
    line = f.readline()
    while line != b"":
        if b"Sum of electronic and thermal Free Energies=" in line:
            gibbs = float(line.split()[-1])
            break
        line = f.readline()
    return float(gibbs)


def SCFOrbitalEnergy(member, tar):
    fh = tar.extractfile(member)

    txt = fh.readlines()
    txt_fwd = tuple([x.strip() for x in txt])
    txt_rev = txt_fwd[::-1]

    occ_energy_levels = list()
    vir_energy_levels = list()

    for i, line in enumerate(txt_rev):
        if line.find(b"Population analysis using the SCF") > -1:
            txt = txt_rev[:i]
            txt = txt[::-1]
            break

    for i, line in enumerate(txt):
        if line.find(b"The electronic state is") > -1:
            txt = txt[i + 1 :]
            break

    for i, line in enumerate(txt):
        if b"Alpha  occ. eigenvalues" in line:
            level = re.findall(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?", line)
            occ_energy_levels.extend(tuple([float(x) for x in level]))
        if b"Alpha virt. eigenvalues" in line:
            level = re.findall(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?", line)
            vir_energy_levels.extend(tuple([float(x) for x in level]))

    occ_energy_levels = np.array(occ_energy_levels)
    vir_energy_levels = np.array(vir_energy_levels)

    homo = occ_energy_levels[-1]
    lumo = vir_energy_levels[0]

    return occ_energy_levels, vir_energy_levels, homo, lumo


# In[18]:


def load_energies(member, tar):
    zpe, e0 = load_zpe_and_scf(member, tar)
    e0_zpe = load_e0_zpe(member, tar)
    gibbs = load_gibbs(member, tar)

    energy = dict()
    energy["scf"] = e0
    energy["zpe_unscaled"] = zpe
    energy["scf_zpe_unscaled"] = e0_zpe
    energy["gibbs"] = gibbs

    return energy


def get_title_card(member, tar, flag=b"Initial command:"):
    f = tar.extractfile(member)
    lines = f.readlines()
    start_inds = [i for i, line in enumerate(lines) if flag in line]
    title_card = b""
    lines = lines[start_inds[-1] :]
    for i, line in enumerate(lines):
        if b" #opt=" in line:
            count = i
            line2 = lines[count]
            while (
                b"------------------------------------------------------" not in line2
            ):
                title_card += line2.strip()
                count += 1
                line2 = lines[count]
            break
    return title_card.decode()


def semiempirical_opt_parser(mol_id, mol_smi, mol_confs_tar=None):

    valid_job = dict()
    failed_job = dict()

    if mol_confs_tar is None:
        ids = mol_id // 1000
        mol_confs_tar = os.path.join(
            "output", "semiempirical_opt", "outputs", f"outputs_{ids}", f"{mol_id}.tar"
        )

    if os.path.isfile(mol_confs_tar):

        valid_job[mol_id] = dict()
        failed_job[mol_id] = dict()

        pre_adj = RDKitMol.FromSmiles(mol_smi).GetAdjacencyMatrix()

        tar = tarfile.open(mol_confs_tar)
        for member in tar:
            conf_id = member.name.split(f"{mol_id}_")[1]
            conf_id = int(conf_id.split(".log")[0])

            job_stat = check_job_status(member, tar)
            if not job_stat:
                failed_job[mol_id][conf_id] = "job status"
                continue

            if not check_freq(member, tar):
                failed_job[mol_id][conf_id] = "freq check"
                continue

            xyz, _, _ = load_geometry(member, tar)
            try:
                post_mol = RDKitMol.FromXYZ(
                    xyz,
                    header=False,
                    sanitize=False,
                )
            except Exception as e:
                failed_job[mol_id][conf_id] = f"rdkit failed with {e}"
                continue
            post_adj = post_mol.GetAdjacencyMatrix()
            if (pre_adj == post_adj).all():

                valid_job[mol_id][conf_id] = dict()
                valid_job[mol_id][conf_id]["mol_smi"] = mol_smi
                valid_job[mol_id][conf_id]["semiempirical_title_card"] = get_title_card(
                    member, tar
                )
                valid_job[mol_id][conf_id]["semiempirical_freq"] = load_freq(
                    member, tar
                )
                (
                    valid_job[mol_id][conf_id]["semiempirical_xyz"],
                    valid_job[mol_id][conf_id]["semiempirical_xyz_dict"],
                    valid_job[mol_id][conf_id]["semiempirical_steps"],
                ) = load_geometry(member, tar)
                (
                    valid_job[mol_id][conf_id]["semiempirical_xyz_std_ori"],
                    valid_job[mol_id][conf_id]["semiempirical_xyz_dict_std_ori"],
                    _,
                ) = load_geometry_std(member, tar)
                valid_job[mol_id][conf_id]["semiempirical_energy"] = load_energies(
                    member, tar
                )
                valid_job[mol_id][conf_id]["semiempirical_cpu"] = get_cpu(member, tar)
                valid_job[mol_id][conf_id]["semiempirical_wall"] = get_wall(member, tar)
            else:
                failed_job[mol_id][conf_id] = "adjacency matrix"
                continue

        if not valid_job[mol_id]:
            del valid_job[mol_id]
            failed_job[mol_id]["reason"] = "all confs failed"
        else:
            if not failed_job[mol_id]:
                del failed_job[mol_id]
    else:
        failed_job[mol_id] = dict()
        failed_job[mol_id]["reason"] = "file not found"

    return failed_job, valid_job


def get_mol_id_to_semiempirical_opted_xyz(valid_jobs):
    mol_id_to_semiempirical_opted_xyz = {}
    for mol_id in valid_jobs:
        ens = np.array(
            [
                conf_dict["semiempirical_energy"]["scf"]
                for conf_id, conf_dict in valid_jobs[mol_id].items()
            ]
        )
        conf_ids = np.array(
            [conf_id for conf_id, conf_dict in valid_jobs[mol_id].items()]
        )
        lowest_conf_ind = conf_ids[np.argsort(ens)[0]]
        xyz = valid_jobs[mol_id][lowest_conf_ind]["semiempirical_xyz_std_ori"]
        xyz = str(len(xyz.splitlines())) + "\n" + f"{mol_id}" + "\n" + xyz
        mol_id_to_semiempirical_opted_xyz[mol_id] = xyz
    return mol_id_to_semiempirical_opted_xyz
