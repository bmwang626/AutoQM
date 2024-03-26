from fileinput import filename
import os
import shutil
import subprocess

from rdkit import Chem
from .file_parser import mol2xyz


def dlpno_sp_calc(mol_id, orca_path, charge, mult, n_procs, job_ram, xyz_DFT_opt):
    sdf = mol_id + '.sdf'
    mol_dir = os.getcwd()

    if xyz_DFT_opt:
        coords = xyz_DFT_opt[mol_id]
    else:
        mol = Chem.SDMolSupplier(sdf, removeHs=False, sanitize=False)[0]
        xyz = mol2xyz(mol)
        coords = "\n".join(xyz.splitlines()[2:])

    #create working directory
    try:
        shutil.rmtree("scratch")
    except:
        pass
    os.makedirs("scratch")
    os.chdir("scratch")

    script = generate_dlpno_sp_input(coords, charge, mult, job_ram, n_procs)

    infile = f"{mol_id}.in"
    with open(infile, "w+") as f:
        f.write(script)

    #run jobs
    orca_command = os.path.join(orca_path, "orca")
    logfile = mol_id + '.log'
    outfile = mol_id + '.out'
    with open(outfile, 'w') as out:
        subprocess.run('{} {} > {}'.format(orca_command, infile, logfile), shell=True, stdout=out, stderr=out)
    
    # check for normal termination
    with open(logfile, "r") as f:
        lines = f.readlines()
    try:
        assert any(["ORCA TERMINATED NORMALLY" in line for line in reversed(lines)])
        os.chdir(mol_dir)
        shutil.copy(os.path.join("scratch", logfile), logfile)
        shutil.rmtree("scratch")
    except:
        os.chdir(mol_dir)
        raise RuntimeError(f"ORCA calculation failed for {mol_id}")

def generate_dlpno_sp_input(level_of_theory: str,
                            xyz_str: str,
                            charge: int,
                            multiplicity: int,
                            memory_mb: int,
                            cpu_threads: int,
                            ) -> str:
    """
    Modified from ACS
    """
    
    script = f"""!{level_of_theory}

%mdci
UseFullLmp2Guess False
end

%maxcore {memory_mb}
%pal
nprocs {cpu_threads}
end

* xyz {charge} {multiplicity}
{xyz_str}
*


"""
    return script

    