#!/usr/bin/env python
# coding: utf-8

import os
import re
import numpy as np

from rdmc.mol import RDKitMol
from rdmc.external.logparser.gaussian import GaussianLog

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


# In[55]:


def read_log_file(self):
    with open(self) as fh:
        txt = fh.readlines()
    log = tuple([x.strip() for x in txt])
    return log


# In[56]:


def split_log(log, flag="Initial command:"):
    indices = [i for i, val in enumerate(log) if val == flag]
    splitted_logs = [log[i:j] for i, j in zip([0] + indices, indices + [None])]
    splitted_logs.pop(0)
    return splitted_logs


# In[57]:


def check_job_status(self):
    for line in reversed(self):
        if "Normal termination" in line:
            return True
        else:
            return False


# In[58]:


def get_cpu(self):
    for line in self:
        if line.find("Job cpu time") > -1:
            days = int(line.split()[3])
            hours = int(line.split()[5])
            mins = int(line.split()[7])
            secs = float(line.split()[9])
            CPU = tuple([days, hours, mins, secs])
            return CPU


# In[59]:


def get_wall(self):
    for line in self:
        if line.find("Elapsed time") > -1:
            days = int(line.split()[2])
            hours = int(line.split()[4])
            mins = int(line.split()[6])
            secs = float(line.split()[8])
            CPU = tuple([days, hours, mins, secs])
            return CPU


# In[85]:


def load_geometry(
    self,
    periodictable=periodictable,
    initial=False,
    input_geom=False,
    standard_orientation=True,
):
    """
    Return the optimum geometry of the molecular configuration from the
    Gaussian log file. If multiple such geometries are identified, only the
    last is returned.
    """
    step = -1
    number, coord, symbol = [], [], []
    if standard_orientation:
        orientation = "Standard orientation:"
    else:
        orientation = "Input orientation:"

    with open(self, "r") as f:
        line = f.readline()
        while line != "":
            # Automatically determine the number of atoms

            if input_geom:
                if "Symbolic Z-matrix:" in line:
                    for i in range(2):
                        line = f.readline()
                    while line.strip() != "":
                        data = line.split()
                        symbol.append(data[0])
                        coord.append([float(data[1]), float(data[2]), float(data[3])])
                        line = f.readline()
                    break

            if orientation in line:
                step += 1
                number, coord = [], []
                for i in range(5):
                    line = f.readline()
                while (
                    "---------------------------------------------------------------------"
                    not in line
                ):
                    data = line.split()
                    number.append(int(data[1]))
                    coord.append([float(data[3]), float(data[4]), float(data[5])])
                    line = f.readline()
            line = f.readline()

            if coord and initial:
                break

    number = np.array(number)
    if not input_geom:
        symbol = [periodictable[x] for x in number]

    xyz_str = make_xyz_str(symbol, coord)
    return xyz_str, step


# In[62]:


def load_freq(self):
    """
    Return the frequencies
    calculation in cm^-1.
    """
    frequencies = []
    with open(self, "r") as f:
        line = f.readline()
        while line != "":
            # Read vibrational frequencies
            if "Frequencies --" in line:
                frequencies.extend(line.split()[2:])
            line = f.readline()

    frequencies = [float(freq) for freq in frequencies]
    frequencies.sort()

    return frequencies


# In[63]:


def check_neg_freq(frequencies, threshold=-50):
    min_freq = min(frequencies)
    has_neg_freq = min_freq < threshold
    if not has_neg_freq:
        min_freq = None
    return has_neg_freq, min_freq


# In[66]:


def load_e0(self):
    with open(self, "r") as f:
        line = f.readline()
        while line != "":
            if "SCF Done:" in line:
                e0 = float(line.split()[4])
            line = f.readline()
    return e0


# In[67]:


def load_zpe(self):
    with open(self, "r") as f:
        s = f.read()
        s = s.replace("\n", "").replace(" ", "")
        zpe = float(re.findall("ZeroPoint=(-*\d+.\d+)", s)[0])
    return zpe


# In[68]:


def load_gibbs(self):
    with open(self, "r") as f:
        line = f.readline()
        while line != "":
            if "Sum of electronic and thermal Free Energies=" in line:
                gibbs = float(line.split()[-1])
                break
            line = f.readline()
    return float(gibbs)


# In[69]:


def process_energy(e0, zpe, zpe_scale_factor):
    e0_zpe = e0 + zpe
    zpe_scaled = zpe * zpe_scale_factor
    e0_zpe_scaled = e0 + zpe_scaled
    return zpe_scaled, zpe_scale_factor, e0_zpe, e0_zpe_scaled


# In[70]:


def load_energies(self, zpe_scale_factor):
    energy = dict()

    e0 = load_e0(self)
    zpe = load_zpe(self)

    energy["scf"] = e0
    energy["zpe_scale_factor"] = zpe_scale_factor
    energy["zpe_unscaled"] = zpe

    composite = process_energy(e0, zpe, zpe_scale_factor)

    energy["zpe_scaled"] = composite[0]
    energy["scf_zpe_unscaled"] = composite[2]
    energy["scf_zpe_scaled"] = composite[3]

    energy["gibbs"] = load_gibbs(self)
    return energy


def dft_opt_freq_parser(
    g16_log,
    is_ts=False,
    check_connectivity=True,
    smi=None,
):
    failed_job = dict()
    valid_job = dict()

    if os.path.isfile(g16_log):
        zpe_scale_factor = 0.986
        # LevelOfTheory(method='wb97xd',basis='def2svp',software='gaussian')": 0.986,  # [4]
        # [4] Calculated as described in 10.1021/ct100326h
        # https://github.com/ReactionMechanismGenerator/RMG-database/blob/main/input/quantum_corrections/data.py

        glog = GaussianLog(g16_log)

        job_stat = glog.success

        if not job_stat:
            failed_job["reason"] = "error termination"
            try:
                failed_job["dft_xyz_std_ori"] = load_geometry(
                    g16_log, standard_orientation=True
                )[0]
                failed_job["dft_initial_xyz_std_ori"] = load_geometry(
                    g16_log, initial=True, standard_orientation=True
                )[0]
                failed_job["dft_xyz_input_ori"] = load_geometry(
                    g16_log, standard_orientation=False
                )[0]
                failed_job["dft_initial_xyz_input_ori"] = load_geometry(
                    g16_log, initial=True, standard_orientation=False
                )[0]
                failed_job["dft_input_xyz"] = load_geometry(g16_log, input_geom=True)[0]
                failed_job["dft_steps"] = load_geometry(g16_log)[1]
                failed_job["dft_cpu"] = get_cpu(read_log_file(g16_log))
                failed_job["dft_wall"] = get_wall(read_log_file(g16_log))
            except:
                failed_job["reason"] = "parser1"
            return failed_job, valid_job

        if check_connectivity:
            pre_adj = RDKitMol.FromSmiles(smi).GetAdjacencyMatrix()

            try:
                post_adj = glog.get_mol(
                    refid=glog.num_all_geoms - 1,  # The last geometry in the job
                    converged=False,
                    sanitize=False,
                    backend="openbabel",
                ).GetAdjacencyMatrix()
            except:
                print(g16_log)
                failed_job["reason"] = "can't get post_adj"
                return failed_job, valid_job

            if not (pre_adj == post_adj).all():
                failed_job["reason"] = "adjacency matrix"
                return failed_job, valid_job

        freqs = load_freq(g16_log)
        has_neg_freq, neg_freq = check_neg_freq(freqs)
        if is_ts:
            pass_freq_check = has_neg_freq
        else:
            pass_freq_check = not has_neg_freq

        if not pass_freq_check:
            failed_job["reason"] = "freq"
            try:
                failed_job["dft_freq"] = freqs
                failed_job["dft_freq_neg"] = has_neg_freq
                failed_job["dft_xyz_std_ori"] = load_geometry(
                    g16_log, standard_orientation=True
                )[0]
                failed_job["dft_initial_xyz_std_ori"] = load_geometry(
                    g16_log, initial=True, standard_orientation=True
                )[0]
                failed_job["dft_xyz_input_ori"] = load_geometry(
                    g16_log, standard_orientation=False
                )[0]
                failed_job["dft_initial_xyz_input_ori"] = load_geometry(
                    g16_log, initial=True, standard_orientation=False
                )[0]
                failed_job["dft_input_xyz"] = load_geometry(g16_log, input_geom=True)[0]
                failed_job["dft_steps"] = load_geometry(g16_log)[1]
                failed_job["dft_cpu"] = get_cpu(read_log_file(g16_log))
                failed_job["dft_wall"] = get_wall(read_log_file(g16_log))
            except:
                failed_job["reason"] = "parser2"

            return failed_job, valid_job

        try:
            valid_job["dft_freq"] = freqs
            valid_job["dft_freq_neg"] = has_neg_freq
            valid_job["dft_xyz_std_ori"] = load_geometry(
                g16_log, standard_orientation=True
            )[0]
            valid_job["dft_initial_xyz_std_ori"] = load_geometry(
                g16_log, initial=True, standard_orientation=True
            )[0]
            valid_job["dft_xyz_input_ori"] = load_geometry(
                g16_log, standard_orientation=False
            )[0]
            valid_job["dft_initial_xyz_input_ori"] = load_geometry(
                g16_log, initial=True, standard_orientation=False
            )[0]
            valid_job["dft_input_xyz"] = load_geometry(g16_log, input_geom=True)[0]
            valid_job["dft_steps"] = load_geometry(g16_log)[1]
            valid_job["dft_cpu"] = get_cpu(read_log_file(g16_log))
            valid_job["dft_wall"] = get_wall(read_log_file(g16_log))
            valid_job["dft_energy"] = load_energies(g16_log, zpe_scale_factor)
        except:
            valid_job = dict()
            failed_job["reason"] = "parser3"
        return failed_job, valid_job
    else:
        failed_job["reason"] = "file not found"
        return failed_job, valid_job
