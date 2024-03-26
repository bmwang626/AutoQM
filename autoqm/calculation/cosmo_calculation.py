from fileinput import filename
import os
import shutil
import subprocess
import csv
import time
import traceback
import pickle as pkl
import tarfile

from .utils import REPLACE_LETTER

from rdkit import Chem
from .file_parser import mol2xyz


def cosmo_calc(
    mol_id,
    cosmotherm_path,
    cosmo_database_path,
    charge,
    mult,
    T_list,
    df_pure,
    xyz,
    scratch_dir,
    tmp_mol_dir,
    save_dir,
    input_dir,
):

    # create and move to working directory
    current_dir = os.getcwd()
    scratch_dir_mol_id = os.path.join(scratch_dir, f"{mol_id}")
    os.makedirs(scratch_dir_mol_id)
    os.chdir(scratch_dir_mol_id)

    # open the tar file
    tar_file = f"{mol_id}.tar"
    tar_file_path = os.path.join(save_dir, tar_file)
    if os.path.exists(tar_file_path):
        tar = tarfile.open(tar_file_path, "r")
        member_basename_list = set(
            os.path.basename(member.name) for member in tar.getmembers()
        )
    else:
        tar = tarfile.open(tar_file_path, "w")
        member_basename_list = set()

    energyfile = f"{mol_id}.energy"
    cosmofile = f"{mol_id}.cosmo"
    if energyfile in member_basename_list and cosmofile in member_basename_list:
        # extract to files
        tar.extract(energyfile, path=".")
        tar.extract(cosmofile, path=".")
        tar.close()
        tar = tarfile.open(tar_file_path, "a")
    else:
        tar.close()
        tar = tarfile.open(tar_file_path, "a")

        # turbomole
        print(f"Running Turbomole for {mol_id}...")

        num_atoms = len(xyz.splitlines())
        xyz = str(num_atoms) + "\n\n" + xyz

        os.makedirs("xyz")
        xyz_mol_id = f"{mol_id}.xyz"
        with open(os.path.join("xyz", xyz_mol_id), "w+") as f:
            f.write(xyz)

        txtfile = f"{mol_id}.txt"
        with open(txtfile, "w+") as f:
            f.write(f"{mol_id} {charge} {mult}")

        # run the job
        logfile = f"{mol_id}.log"
        outfile = f"{mol_id}.out"
        with open(outfile, "w") as out:
            subprocess.run(
                f"calculate -l {txtfile} -m BP-TZVPD-FINE-COSMO-SP -f xyz -din xyz > {logfile}",
                shell=True,
                stdout=out,
                stderr=out,
            )
            subprocess.run(
                f"calculate -l {txtfile} -m BP-TZVPD-GAS-SP -f xyz -din xyz > {logfile}",
                shell=True,
                stdout=out,
                stderr=out,
            )

        cosmo_done = False
        energy_done = False

        # copy the cosmo and energy files
        for file in os.listdir("CosmofilesBP-TZVPD-FINE-COSMO-SP"):
            if file.endswith("cosmo"):
                shutil.copyfile(
                    os.path.join("CosmofilesBP-TZVPD-FINE-COSMO-SP", file), file
                )
                tar.add(file)
                tar.close()
                tar = tarfile.open(tar_file_path, "a")
                cosmo_done = True
                break

        for file in os.listdir("EnergyfilesBP-TZVPD-FINE-COSMO-SP"):
            if file.endswith("energy"):
                shutil.copyfile(
                    os.path.join("EnergyfilesBP-TZVPD-FINE-COSMO-SP", file), file
                )
                tar.add(file)
                tar.close()
                tar = tarfile.open(tar_file_path, "a")
                energy_done = True
                break

        if not (cosmo_done and energy_done):
            shutil.copyfile(
                os.path.join("xyz", xyz_mol_id), os.path.join(tmp_mol_dir, xyz_mol_id)
            )
            shutil.copyfile(txtfile, os.path.join(tmp_mol_dir, txtfile))
            shutil.copyfile(outfile, os.path.join(tmp_mol_dir, outfile))
            shutil.copyfile(logfile, os.path.join(tmp_mol_dir, logfile))
            print(f"Turbomole calculation failed for {mol_id}")
            print("Output files:")
            with open(outfile, "r") as f:
                print(f.read())
            print("Log files:")
            with open(logfile, "r") as f:
                print(f.read())
            os.chdir(current_dir)
            return

        print(f"Turbomole calculation done for {mol_id}")

    # prepare for cosmo calculation
    for index, row in df_pure.iterrows():
        cosmo_name = "".join(
            letter if letter not in REPLACE_LETTER else REPLACE_LETTER[letter]
            for letter in row.cosmo_name
        )
        inpfile = f"{mol_id}_{cosmo_name}.inp"
        tabfile = f"{mol_id}_{cosmo_name}.tab"
        outfile = f"{mol_id}_{cosmo_name}.out"

        if tabfile in member_basename_list:
            continue

        print(f"Running COSMO calculation for {mol_id} in {index} {row.cosmo_name}...")

        script = generate_cosmo_input(
            str(mol_id), cosmotherm_path, cosmo_database_path, T_list, row
        )

        with open(inpfile, "w+") as f:
            f.write(script)

        cosmo_command = os.path.join(
            cosmotherm_path, "COSMOtherm", "BIN-LINUX", "cosmotherm"
        )
        subprocess.run(f"{cosmo_command} {inpfile}", shell=True)

        if not os.path.exists(tabfile):
            shutil.copyfile(inpfile, os.path.join(tmp_mol_dir, inpfile))
            shutil.copyfile(outfile, os.path.join(tmp_mol_dir, outfile))
            print(f"COSMO calculation failed for {mol_id} in {index} {row.cosmo_name}")
            with open(inpfile, "r") as f:
                print(f.read())
            with open(outfile, "r") as f:
                print(f.read())
            os.chdir(current_dir)
            return
        else:
            tar.add(inpfile)
            tar.add(tabfile)
            tar.add(outfile)
            tar.close()
            tar = tarfile.open(tar_file_path, "a")

            print(f"COSMO calculation done for {mol_id} in {index} {row.cosmo_name}")

    tar.close()

    print("Removing temporary folder...")
    if os.path.exists(tmp_mol_dir):
        if not os.listdir(tmp_mol_dir):
            shutil.rmtree(tmp_mol_dir)
        else:
            print(f"Some jobs for {mol_id} failed. See {tmp_mol_dir} for details.")
    else:
        print("Removed by other worker? Skipping...")

    print("Removing temporary file...")
    try:
        os.remove(os.path.join(input_dir, f"{mol_id}.tmp"))
    except FileNotFoundError as e:
        print(e)
        print("Removed by other worker? Skipping...")

    os.chdir(current_dir)
    shutil.rmtree(scratch_dir_mol_id)


def generate_cosmo_input(name, cosmotherm_path, cosmo_database_path, T_list, row):
    """
    Modified from ACS and Yunsie's code
    """

    script = f"""ctd = BP_TZVPD_FINE_21.ctd cdir = "{cosmotherm_path}/COSMOthermX/../COSMOtherm/CTDATA-FILES" ldir = "{cosmotherm_path}/COSMOthermX/../licensefiles"
notempty wtln ehfile
!! generated by COSMOthermX !!
"""

    # solvent
    first_letter = row.cosmo_name[0]
    if not first_letter.isalpha() and not first_letter.isnumeric():
        first_letter = "0"
    if row.source == "COSMOtherm":
        solvent_dir = (
            f"{cosmotherm_path}/COSMOtherm/DATABASE-COSMO/BP-TZVPD-FINE/{first_letter}"
        )
    elif row.source == "COSMObase":
        solvent_dir = f"{cosmo_database_path}/BP-TZVPD-FINE/{first_letter}"
    script += 'f = "' + row.cosmo_name + '_c0.cosmo" fdir="' + solvent_dir + '"'
    if int(row.cosmo_conf) > 1:
        script += ' Comp = "' + row.cosmo_name + '" [ VPfile'
        for k in range(1, int(row.cosmo_conf)):
            script += (
                '\nf = "'
                + row.cosmo_name
                + "_c"
                + str(k)
                + '.cosmo" fdir="'
                + solvent_dir
                + '"'
            )
        script += " ]\n"
    else:
        script += " VPfile\n"

    # solute
    script += 'f = "' + name + '.cosmo" fdir="." VPfile\n'
    for T in T_list:
        script += "henry  xh={ 1 0 } tk=" + str(T) + " GSOLV  \n"
    return script


def read_cosmo_tab_result(tab_file_path):
    """
    Modified from Yunsie's code
    """
    each_data_list = []
    # initialize everything
    solvent_name, solute_name, temp = None, None, None
    result_values = None
    with open(tab_file_path, "r") as f:
        line = f.readline()
        while line:
            # get the temperature and mole fraction
            if "Settings  job" in line:
                temp = line.split("T=")[1].split("K")[0].strip()  # temp in K

            # get the result values
            if "Nr Compound" in line:
                line = f.readline()
                solvent_name = line.split()[1]
                line = f.readline()
                solute_name = line.split()[1]
                result_values = line.split()[
                    2:6
                ]  # H (in bar), ln(gamma), pv (vapor pressure in bar), Gsolv (kcal/mol)
                # save the result as one list
                each_data_list.append(
                    [solvent_name, solute_name, temp] + result_values + [None]
                )
                # initialize everything
                solvent_name, solute_name, temp = None, None, None
                result_values = None
            line = f.readline()
    return each_data_list


def get_dHsolv_value(each_data_list):
    # compute solvation enthalpy
    dGsolv_temp_dict = {}
    ind_298 = None
    for z in range(len(each_data_list)):
        temp = each_data_list[z][2]
        dGsolv = each_data_list[z][6]
        dGsolv_temp_dict[temp] = dGsolv
        if temp == "298.15":
            ind_298 = z
    dGsolv_298 = float(dGsolv_temp_dict["298.15"])
    dSsolv_298 = -(
        float(dGsolv_temp_dict["299.15"]) - float(dGsolv_temp_dict["297.15"])
    ) / (299.15 - 297.15)
    dHsolv_298 = dGsolv_298 + 298.15 * dSsolv_298
    each_data_list[ind_298][7] = "%.8f" % dHsolv_298
    return each_data_list
