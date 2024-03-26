from genericpath import isfile
from rdkit import Chem
import os
import shutil
import subprocess
import traceback
import tarfile

import numpy as np

from .log_parser import XtbLog, G16Log
from .file_parser import mol2xyz, xyz2com, write_mol_to_sdf, write_mols_to_sdf


def run_xtb_opt(
    xyz, charge, mult, mol_id, rdmc_path, g16_path, n_procs, job_ram, level_of_theory
):
    comfile = f"{mol_id}.gjf"
    logfile = f"{mol_id}.log"
    outfile = f"{mol_id}.out"

    g16_command = os.path.join(g16_path, "g16")

    head = '%nprocshared={}\n%mem={}mb\n{}\nexternal="{}/rdmc/external/xtb_tools/xtb_gaussian.pl --gfn 2 -P"\n'.format(
        n_procs, job_ram, level_of_theory, rdmc_path
    )

    xyz2com(xyz, head=head, comfile=comfile, charge=charge, mult=mult, footer="\n")

    with open(outfile, "w") as out:
        subprocess.run(
            "{} < {} >> {}".format(g16_command, comfile, logfile),
            shell=True,
            stdout=out,
            stderr=out,
        )


def semiempirical_opt(
    mol_id,
    charge,
    mult,
    xyz_FF_dict,
    xtb_path,
    rdmc_path,
    g16_path,
    level_of_theory,
    n_procs,
    job_ram,
    scratch_dir,
    tmp_mol_dir,
    suboutputs_dir,
    subinputs_dir,
):
    current_dir = os.getcwd()

    for conf_ind, xyz in xyz_FF_dict[mol_id].items():
        comfile = f"{mol_id}_{conf_ind}.gjf"
        logfile = f"{mol_id}_{conf_ind}.log"
        outfile = f"{mol_id}_{conf_ind}.out"

        if os.path.exists(os.path.join(tmp_mol_dir, logfile)):
            continue

        conf_scratch_dir = os.path.join(scratch_dir, f"{mol_id}_{conf_ind}")
        os.makedirs(conf_scratch_dir)
        os.chdir(conf_scratch_dir)

        run_xtb_opt(
            xyz,
            charge,
            mult,
            f"{mol_id}_{conf_ind}",
            rdmc_path,
            g16_path,
            n_procs,
            job_ram,
            level_of_theory,
        )
        shutil.copyfile(logfile, os.path.join(tmp_mol_dir, logfile))
        os.chdir(current_dir)

    mol_scratch_dir = os.path.join(scratch_dir, f"{mol_id}")
    os.makedirs(mol_scratch_dir)
    os.chdir(mol_scratch_dir)

    # tar the log files
    tar_file = f"{mol_id}.tar"
    tar = tarfile.open(tar_file, "w")
    for conf_ind, xyz in xyz_FF_dict[mol_id].items():
        logfile = f"{mol_id}_{conf_ind}.log"
        tar.add(os.path.join(tmp_mol_dir, logfile))
    tar.close()

    shutil.copy(tar_file, os.path.join(suboutputs_dir, tar_file))
    try:
        os.remove(os.path.join(subinputs_dir, f"{mol_id}.tmp"))
    except FileNotFoundError as e:
        print(e)
        pass
    shutil.rmtree(tmp_mol_dir)
    os.chdir(current_dir)


def xtb_status(folder, molid):

    try:
        log = XtbLog(os.path.join(folder, "{}_freq.log".format(molid)))
    except:
        raise RuntimeError(f"xtb log file not found for {molid}")

    if log.termination:
        peaks = log.wavenum
        if np.min(peaks) < 0:
            raise RuntimeError("imaginary frequency found for {}".format(molid))
        else:
            return "{}_opt.sdf".format(molid)
    else:
        raise RuntimeError("xtb optimization did not finish for {}".format(molid))
