import os
import tarfile
import shutil
import rdkit.Chem as Chem

from rdmc.mol import RDKitMol
from rdmc.forcefield import OpenBabelFF
from rdmc.ts import get_formed_and_broken_bonds

from autoqm.calculation.semiempirical_calculation import run_xtb_opt
from autoqm.calculation.utils import mol2charge, mol2mult, mol2xyz


def reset_r_p_complex_semi_opt(
    rxn_smi,
    ts_xyz,
    ts_id,
    rdmc_path,
    g16_path,
    level_of_theory,
    n_procs,
    job_ram,
    subinputs_dir,
    suboutputs_dir,
    scratch_dir,
):
    current_dir = os.path.abspath(os.getcwd())

    r_complex_smi, p_complex_smi = rxn_smi.split(">>")
    r_complex = RDKitMol.FromSmiles(r_complex_smi)
    p_complex = RDKitMol.FromSmiles(p_complex_smi)
    ts_mol = RDKitMol.FromXYZ(ts_xyz, header=False, sanitize=False)

    formed_bonds, broken_bonds = get_formed_and_broken_bonds(r_complex, p_complex)

    r_complex_id = f"{ts_id}_r"
    rmol_scratch_dir = os.path.join(scratch_dir, r_complex_id)
    os.makedirs(rmol_scratch_dir)
    os.chdir(rmol_scratch_dir)
    new_r_complex = reset_r_complex(ts_mol, r_complex, formed_bonds)
    xyz = mol2xyz(new_r_complex)
    charge = mol2charge(new_r_complex)
    mult = mol2mult(new_r_complex)
    run_xtb_opt(
        xyz,
        charge,
        mult,
        r_complex_id,
        rdmc_path,
        g16_path,
        n_procs,
        job_ram,
        level_of_theory,
    )
    # shutil.copyfile(f"{r_complex_id}.gjf", os.path.join(suboutputs_dir, f"{r_complex_id}.gjf"))
    # shutil.copyfile(f"{r_complex_id}.log", os.path.join(suboutputs_dir, f"{r_complex_id}.log"))
    # shutil.copyfile(f"{r_complex_id}.out", os.path.join(suboutputs_dir, f"{r_complex_id}.out"))
    os.chdir(current_dir)

    p_complex_id = f"{ts_id}_p"
    pmol_scratch_dir = os.path.join(scratch_dir, p_complex_id)
    os.makedirs(pmol_scratch_dir)
    os.chdir(pmol_scratch_dir)
    new_p_complex = reset_p_complex(new_r_complex, p_complex, broken_bonds)
    xyz = mol2xyz(new_p_complex)
    charge = mol2charge(new_p_complex)
    mult = mol2mult(new_p_complex)
    run_xtb_opt(
        xyz,
        charge,
        mult,
        p_complex_id,
        rdmc_path,
        g16_path,
        n_procs,
        job_ram,
        level_of_theory,
    )
    # shutil.copyfile(f"{p_complex_id}.gjf", os.path.join(suboutputs_dir, f"{p_complex_id}.gjf"))
    # shutil.copyfile(f"{p_complex_id}.log", os.path.join(suboutputs_dir, f"{p_complex_id}.log"))
    # shutil.copyfile(f"{p_complex_id}.out", os.path.join(suboutputs_dir, f"{p_complex_id}.out"))
    os.chdir(current_dir)

    ts_scratch_dir = os.path.join(scratch_dir, ts_id)
    os.makedirs(ts_scratch_dir)
    os.chdir(ts_scratch_dir)

    # tar the cosmo, energy and tab files
    tar_file = f"{ts_id}.tar"
    tar = tarfile.open(tar_file, "w")
    tar.add(os.path.join(rmol_scratch_dir, f"{r_complex_id}.log"))
    tar.add(os.path.join(pmol_scratch_dir, f"{p_complex_id}.log"))
    tar.close()

    shutil.copyfile(tar_file, os.path.join(suboutputs_dir, tar_file))
    os.chdir(current_dir)
    try:
        os.remove(os.path.join(subinputs_dir, f"{ts_id}.tmp"))
    except FileNotFoundError:
        print("File not found")
        print(os.path.join(subinputs_dir, f"{ts_id}.tmp"))
    shutil.rmtree(rmol_scratch_dir)
    shutil.rmtree(pmol_scratch_dir)
    shutil.rmtree(ts_scratch_dir)


def reset_r_p_complex_ff_opt(
    rxn_smi, ts_xyz, ts_id, subinputs_dir, suboutputs_dir, scratch_dir
):
    current_dir = os.path.abspath(os.getcwd())

    r_complex_smi, p_complex_smi = rxn_smi.split(">>")
    r_complex = RDKitMol.FromSmiles(r_complex_smi, removeHs=False, sanitize=False)
    p_complex = RDKitMol.FromSmiles(p_complex_smi, removeHs=False, sanitize=False)
    ts_mol = RDKitMol.FromXYZ(ts_xyz, header=False, sanitize=False)

    formed_bonds, broken_bonds = get_formed_and_broken_bonds(r_complex, p_complex)

    # get bond length of reacting bonds
    ts_conf = ts_mol.GetConformer()
    reacting_bond_lengths = [
        ts_conf.GetBondLength(b) for b in broken_bonds + formed_bonds
    ]
    if any([b > 2.0 for b in reacting_bond_lengths]):
        print(
            f"Warning: reacting bond length of {ts_id} {rxn_smi} is greater than 2.0 Angstrom: {reacting_bond_lengths}. Skip this TS."
        )
        try:
            os.remove(os.path.join(subinputs_dir, f"rxn_{ts_id}.tmp"))
        except FileNotFoundError:
            print("File not found. Continuing...")
            print(os.path.join(subinputs_dir, f"rxn_{ts_id}.tmp"))
        return

    r_complex_id = f"rxn_{ts_id}_r"
    rmol_scratch_dir = os.path.join(scratch_dir, r_complex_id)
    os.makedirs(rmol_scratch_dir)
    os.chdir(rmol_scratch_dir)
    new_r_complex = reset_r_complex(ts_mol, r_complex, formed_bonds)
    os.chdir(current_dir)

    p_complex_id = f"rxn_{ts_id}_p"
    pmol_scratch_dir = os.path.join(scratch_dir, p_complex_id)
    os.makedirs(pmol_scratch_dir)
    os.chdir(pmol_scratch_dir)
    new_p_complex = reset_p_complex(new_r_complex, p_complex, broken_bonds)

    new_r_complex._mol.SetProp("_Name", r_complex_smi)
    new_p_complex._mol.SetProp("_Name", p_complex_smi)
    ts_mol._mol.SetProp("_Name", rxn_smi)

    sdf_file = f"rxn_{ts_id}.sdf"
    writer = Chem.rdmolfiles.SDWriter(sdf_file)
    writer.write(new_r_complex._mol)
    writer.write(new_p_complex._mol)
    writer.write(ts_mol._mol)
    writer.close()

    shutil.copyfile(sdf_file, os.path.join(suboutputs_dir, sdf_file))
    os.chdir(current_dir)

    try:
        os.remove(os.path.join(subinputs_dir, f"rxn_{ts_id}.tmp"))
    except FileNotFoundError:
        print("File not found. Continuing...")
        print(os.path.join(subinputs_dir, f"rxn_{ts_id}.tmp"))
    shutil.rmtree(rmol_scratch_dir)
    shutil.rmtree(pmol_scratch_dir)


def reset_r_complex(ts_mol, r_complex, formed_bonds):
    # copy current r_complex and set new positions
    new_r_complex = r_complex.Copy(quickCopy=True)
    new_r_complex.SetPositions(ts_mol.GetPositions())

    # setup first minimization with broken bond constraints
    obff = OpenBabelFF(force_field="uff")
    obff.setup(new_r_complex)
    ts_conf = ts_mol.GetConformer()
    current_distances = [ts_conf.GetBondLength(b) for b in formed_bonds]
    for b, d in zip(formed_bonds, current_distances):
        obff.add_distance_constraint(b, 1.5 * d)
    obff.optimize(max_step=2000)

    # second minimization without constraints
    obff.constraints = None
    obff.optimize(max_step=2000)
    new_r_complex = obff.get_optimized_mol()

    # third optimization with MMFF94s
    obff = OpenBabelFF(force_field="mmff94s")
    obff.setup(new_r_complex)
    obff.optimize(max_step=2000)
    new_r_complex = obff.get_optimized_mol()

    return new_r_complex


def reset_p_complex(new_r_complex, p_complex, broken_bonds):
    # copy current p_complex and set new positions
    new_p_complex = p_complex.Copy(quickCopy=True)
    new_p_complex.SetPositions(new_r_complex.GetPositions())

    # setup first minimization with broken bond constraints
    obff = OpenBabelFF(force_field="uff")
    obff.setup(new_p_complex)
    r_conf = new_r_complex.GetConformer()
    current_distances = [r_conf.GetBondLength(b) for b in broken_bonds]
    for b, d in zip(broken_bonds, current_distances):
        obff.add_distance_constraint(b, 1.5 * d)
    obff.optimize(max_step=2000)

    # second minimization without constraints
    obff.constraints = None
    obff.optimize(max_step=2000)
    new_p_complex = obff.get_optimized_mol()

    # third optimization with MMFF94s
    obff = OpenBabelFF(force_field="mmff94s")
    obff.setup(new_p_complex)
    obff.optimize(max_step=2000)
    new_p_complex = obff.get_optimized_mol()

    return new_p_complex
