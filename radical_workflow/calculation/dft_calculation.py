import shutil
from rdkit import Chem
import copy
import csv
import os
import subprocess
import numpy as np
import time

from .file_parser import mol2xyz, xyz2com, write_mol_to_sdf
from .grab_QM_descriptors import read_log
from .log_parser import G16Log
from radical_workflow.parser.dft_opt_freq_parser import dft_opt_freq_parser


def dft_scf_qm_descriptor(
    folder, sdf, g16_path, level_of_theory, n_procs, logger, job_ram, base_charge
):
    basename = os.path.basename(sdf)

    parent_folder = os.getcwd()
    os.chdir(folder)

    try:
        mol_id = os.path.splitext(basename)[0]

        xyz = mol2xyz(Chem.SDMolSupplier(sdf, removeHs=False, sanitize=False)[0])

        pwd = os.getcwd()

        g16_command = os.path.join(g16_path, "g16")
        QM_descriptors = {}
        for jobtype in ["neutral", "plus1", "minus1"]:
            os.makedirs(jobtype, exist_ok=True)

            if jobtype == "neutral":
                charge = base_charge
                mult = 1
                head = (
                    "%chk={}.chk\n%nprocshared={}\n%mem={}mb\n# b3lyp/def2svp nmr=GIAO scf=(maxcycle=512, xqc) "
                    "pop=(full,mbs,hirshfeld,nbo6read)\n".format(
                        mol_id, n_procs, job_ram
                    )
                )
            elif jobtype == "plus1":
                charge = base_charge + 1
                mult = 2
                head = (
                    "%chk={}.chk\n%nprocshared={}\n%mem={}mb\n# b3lyp/def2svp scf=(maxcycle=512, xqc) "
                    "pop=(full,mbs,hirshfeld,nbo6read)\n".format(
                        mol_id, n_procs, job_ram
                    )
                )
            elif jobtype == "minus1":
                charge = base_charge - 1
                mult = 2
                head = (
                    "%chk={}.chk\n%nprocshared={}\n%mem={}mb\n# b3lyp/def2svp scf=(maxcycle=512, xqc) "
                    "pop=(full,mbs,hirshfeld,nbo6read)\n".format(
                        mol_id, n_procs, job_ram
                    )
                )

            os.chdir(jobtype)
            comfile = mol_id + ".gjf"
            xyz2com(
                xyz,
                head=head,
                comfile=comfile,
                charge=charge,
                mult=mult,
                footer="$NBO BNDIDX $END\n",
            )

            logfile = mol_id + ".log"
            outfile = mol_id + ".out"
            if not os.path.exists(outfile):
                with open(outfile, "w") as out:
                    subprocess.run(
                        "{} < {} >> {}".format(g16_command, comfile, logfile),
                        shell=True,
                        stdout=out,
                        stderr=out,
                    )
                    QM_descriptors[jobtype] = read_log(logfile, jobtype)
            else:
                with open(outfile) as f:
                    if "Aborted" in f.read():
                        with open(outfile, "w") as out:
                            subprocess.run(
                                "{} < {} >> {}".format(g16_command, comfile, logfile),
                                shell=True,
                                stdout=out,
                                stderr=out,
                            )
                            QM_descriptors[jobtype] = read_log(logfile, jobtype)
                    else:
                        QM_descriptors[jobtype] = read_log(logfile, jobtype)
            os.chdir(pwd)

        QM_descriptors_return = copy.deepcopy(QM_descriptors)
        QM_descriptor_calc = dict()

        # charges and fukui indices
        for charge in ["mulliken_charge", "hirshfeld_charges", "NPA_Charge"]:
            QM_descriptor_calc["{}_plus1".format(charge)] = QM_descriptors["plus1"][
                charge
            ]
            QM_descriptor_calc["{}_minus1".format(charge)] = QM_descriptors["minus1"][
                charge
            ]

            QM_descriptor_calc["{}_fukui_elec".format(charge)] = (
                QM_descriptors["neutral"][charge] - QM_descriptors["minus1"][charge]
            )
            QM_descriptor_calc["{}_fukui_neu".format(charge)] = (
                QM_descriptors["plus1"][charge] - QM_descriptors["neutral"][charge]
            )

        # spin density
        for spin in ["mulliken_spin_density", "hirshfeld_spin_density"]:
            QM_descriptor_calc["{}_plus1".format(spin)] = QM_descriptors["plus1"][spin]
            QM_descriptor_calc["{}_minus1".format(charge)] = QM_descriptors["minus1"][
                spin
            ]

        # SCF
        QM_descriptor_calc["SCF_plus1"] = QM_descriptors["plus1"]["SCF"]
        QM_descriptor_calc["SCF_minus1"] = QM_descriptors["minus1"]["SCF"]

        QM_descriptors_return["calculated"] = copy.deepcopy(QM_descriptor_calc)

        os.remove(sdf)
    finally:
        os.chdir(parent_folder)

    return QM_descriptors_return


def dft_scf_opt(
    mol_id,
    xyz,
    g16_path,
    level_of_theory,
    n_procs,
    job_ram,
    base_charge,
    mult,
    scratch_dir,
    input_dir,
    output_dir,
):
    current_dir = os.getcwd()

    mol_scratch_dir = os.path.join(scratch_dir, f"{mol_id}")
    print(mol_scratch_dir)

    os.makedirs(mol_scratch_dir)
    os.chdir(mol_scratch_dir)

    g16_command = os.path.join(g16_path, "g16")
    head = "%chk={}.chk\n%nprocshared={}\n%mem={}mb\n{}\n".format(
        mol_id, n_procs, job_ram, level_of_theory
    )

    comfile = f"{mol_id}.gjf"
    xyz2com(xyz, head=head, comfile=comfile, charge=base_charge, mult=mult, footer="\n")
    shutil.copyfile(comfile, os.path.join(output_dir, f"{mol_id}.gjf"))

    logfile = f"{mol_id}.log"
    outfile = f"{mol_id}.out"

    start_time = time.time()
    with open(outfile, "w") as out:
        subprocess.run(
            "{} < {} >> {}".format(g16_command, comfile, logfile),
            shell=True,
            stdout=out,
            stderr=out,
        )
    end_time = time.time()
    print(
        f"Optimization of {mol_id} with {level_of_theory} took {end_time - start_time} seconds."
    )

    shutil.copyfile(outfile, os.path.join(output_dir, f"{mol_id}.out"))
    shutil.copyfile(logfile, os.path.join(output_dir, f"{mol_id}.log"))
    os.remove(os.path.join(input_dir, f"{mol_id}.in"))

    os.chdir(current_dir)


def dft_scf_sp(
    mol_id, g16_path, level_of_theory, n_procs, logger, job_ram, base_charge, mult
):
    sdf = mol_id + ".sdf"

    mol = Chem.SDMolSupplier(sdf, removeHs=False, sanitize=False)[0]
    xyz = mol2xyz(mol)

    g16_command = os.path.join(g16_path, "g16")
    head = "%chk={}.chk\n%nprocshared={}\n%mem={}mb\n{}\n".format(
        mol_id, n_procs, job_ram, level_of_theory
    )

    comfile = mol_id + ".gjf"
    xyz2com(xyz, head=head, comfile=comfile, charge=base_charge, mult=mult, footer="\n")

    logfile = mol_id + ".log"
    outfile = mol_id + ".out"
    with open(outfile, "w") as out:
        subprocess.run(
            "{} < {} >> {}".format(g16_command, comfile, logfile),
            shell=True,
            stdout=out,
            stderr=out,
        )

    os.remove(sdf)


def save_dft_sp_results(
    folder, done_jobs_record, task_id, mol_id_to_smi_dict, semiempirical_methods
):
    """extract dft single point calculation results for geometry optimized with different semiempirical methods"""
    result_file_path = f"dft_sp_result_{task_id}.csv"
    header = ["id", "smiles", "GFN2-XTB_opt_DFT_sp", "am1_opt_DFT_sp", "pm7_opt_DFT_sp"]
    with open(result_file_path, "w") as csvfile:
        # creating a csv writer object
        csvwriter = csv.writer(csvfile)
        # writing the header
        csvwriter.writerow(header)

        for mol_id, semiempirical_methods in done_jobs_record.test_DFT_sp.items():
            each_data_list = [mol_id, mol_id_to_smi_dict[mol_id]]
            for semiempirical_method in semiempirical_methods:
                log_file_path = os.path.join(
                    folder, mol_id, semiempirical_method, mol_id + ".log"
                )
                g16log = G16Log(log_file_path)
                en = g16log.E
                each_data_list.append(str(en))
            csvwriter.writerow(each_data_list)
