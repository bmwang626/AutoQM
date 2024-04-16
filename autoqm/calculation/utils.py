from pathlib import Path

REPLACE_LETTER = {"(": "_", ")": "_", "'": "_"}


def mol2xyz(mol):
    return mol.ToXYZ()


def mol2charge(mol):
    return mol.GetFormalCharge()


def mol2mult(mol):
    num_radical_elec = 0
    for atom in mol.GetAtoms():
        num_radical_elec += atom.GetNumRadicalElectrons()
    return num_radical_elec + 1


def add_shared_arguments(parser):
    input_parser = parser.add_argument_group("Input")
    input_parser.add_argument(
        "--input_file",
        required=True,
        type=Path,
        help="input CSV file containing the species information",
    )
    input_parser.add_argument(
        "--id_column", required=True, help="column name for the species id (needs to be integer)."
    )
    input_parser.add_argument(
        "--smiles_column", required=True, help="column name for the SMILES string"
    )
    input_parser.add_argument(
        "--xyz_column", required=True, help="column name for the xyz string"
    )
    input_parser.add_argument(
        "--xyz_DFT_opt_dict",
        default=None,
        type=Path,
        help="pickle file containing a dictionary to map between the mol_id and DFT-optimized xyz for following calculations",
    )

    parser.add_argument(
        "--scratch_dir", required=True, type=Path, help="scfratch directory"
    )
    parser.add_argument(
        "--task_id",
        type=int,
        default=0,
        help="task id for the calculation",
    )
    parser.add_argument(
        "--num_tasks",
        type=int,
        default=1,
        help="number of tasks for the calculation",
    )
    return parser


def add_cosmo_arguments(parser):
    parser.add_argument(
        "--COSMO_temperatures",
        nargs="+",
        required=False,
        default=["297.15", "298.15", "299.15"],
        help="temperatures used for COSMO calculation",
    )
    parser.add_argument(
        "--COSMO_input_pure_solvents",
        required=True,
        help="input file containing pure solvents used for COSMO calculation.",
    )
    parser.add_argument("--COSMOtherm_path", required=True, help="path to COSMOthermo")
    parser.add_argument(
        "--COSMO_database_path", required=True, type=Path, help="path to COSMO_database"
    )


def add_xtb_arguments(parser):
    parser.add_argument(
        "--XTB_path", required=True, type=Path, help="path to installed XTB"
    )
    parser.add_argument(
        "--RDMC_path",
        required=True,
        type=Path,
        help="path to RDMC to use xtb-gaussian script for xtb optimization calculation.",
    )
    return parser


def add_orca_arguments(parser):
    parser.add_argument("--ORCA_path", required=True, type=Path, help="path to ORCA")
    return parser


def add_gaussian_arguments(parser):
    gaussian_parser = parser.add_argument_group("Gaussian arguments")
    gaussian_parser.add_argument(
        "--template_file",
        required=True,
        help="template file for Gaussian input file",
    )
    gaussian_parser.add_argument(
        "--g16_path", required=True, type=Path, help="path to installed Gaussian 16"
    )
    return parser
