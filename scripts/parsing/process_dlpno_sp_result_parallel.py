import os
import sys
import pickle as pkl
import pandas as pd
from joblib import Parallel, delayed


class OrcaLog(object):
    def __init__(self, path):
        self.path = path

    def check_for_errors(self):
        """
        Checks for common errors in an Orca log file.
        If any are found, this method will raise an error and crash.
        """
        with open(os.path.join(self.path), "r") as f:
            lines = f.readlines()
            error = None
            for line in reversed(lines):
                # check for common error messages
                if "This wavefunction IS NOT CONVERGED!" in line:
                    error = "This wavefunction IS NOT CONVERGED!"
                    break
                elif "ORCA finished by error termination in SCF" in line:
                    error = "SCF"
                    break
                elif "ORCA finished by error termination in MDCI" in line:
                    error = "MDCI"
                    break
                elif "Error : multiplicity" in line:
                    error = "The multiplicity and charge combination are wrong."
                    break
                elif "Wavefunction not fully converged!" in line:
                    error = "Wavefunction not fully converged!"
                    break
                elif "ORCA finished by error termination in GTOInt" in line:
                    error = "GTOInt"
                    break
                elif "ORCA TERMINATED NORMALLY" in line:
                    break
        return error

    def load_energy(self):
        """
        Load the energy in J/ml from an Orca log file. Only the last energy
        in the file is returned. The zero-point energy is *not* included in
        the returned value.
        """
        e_elect = None
        with open(self.path, "r") as f:
            for line in f:
                if "FINAL SINGLE POINT ENERGY" in line:  # for all methods in Orca
                    try:
                        e_elect = float(line.split()[-1])
                    except:
                        print(self.path)
                        print(line)
                        raise
        if e_elect is None:
            raise LogError("Unable to find energy in Orca output file.")
        return e_elect


def parser(mol_id, mol_smi):

    ids = str(int(int(mol_id.split("id")[1]) / 1000))
    orca_log = os.path.join(
        "output", "DLPNO_sp", "outputs", f"outputs_{ids}", f"{mol_id}.log"
    )
    failed_jobs = dict()
    valid_job = dict()

    if os.path.isfile(orca_log):

        olog = OrcaLog(orca_log)

        error = olog.check_for_errors()

        if error is not None:
            failed_jobs[mol_id] = dict()
            failed_jobs[mol_id]["status"] = False
            failed_jobs[mol_id]["reason"] = error
            return failed_jobs, valid_job

        valid_job[mol_id] = dict()
        valid_job[mol_id]["mol_smi"] = mol_smi
        valid_job[mol_id]["dlpno_energy"] = olog.load_energy()

    else:
        failed_jobs[mol_id] = dict()
        failed_jobs[mol_id]["reason"] = "file not found"
    return failed_jobs, valid_job


def main(input_smiles_path, output_file_name, n_jobs):

    df = pd.read_csv(input_smiles_path)
    mol_ids = df["id"].tolist()
    mol_id_to_smi = dict(zip(df["id"].tolist(), df["smiles"].tolist()))

    out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(
        delayed(parser)(mol_id, mol_id_to_smi[mol_id]) for mol_id in mol_ids
    )

    failed_jobs = dict()
    valid_jobs = dict()
    for failed_job, valid_job in out:
        failed_jobs.update(failed_job)
        valid_jobs.update(valid_job)

    with open(os.path.join(f"{output_file_name}.pkl"), "wb") as outfile:
        pkl.dump(valid_jobs, outfile, protocol=pkl.HIGHEST_PROTOCOL)

    with open(os.path.join(f"{output_file_name}_failed.pkl"), "wb") as outfile:
        pkl.dump(failed_jobs, outfile, protocol=pkl.HIGHEST_PROTOCOL)

    print(f"Total number of jobs: {len(mol_ids)}")
    print(f"Number of failed jobs: {len(failed_jobs)}")
    print(failed_jobs)


if __name__ == "__main__":

    input_smiles_path = sys.argv[1]
    output_file_name = sys.argv[2]
    n_jobs = int(sys.argv[3])

    main(input_smiles_path, output_file_name, n_jobs)
