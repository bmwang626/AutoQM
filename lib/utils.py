import logging
import os
import yaml
from rdkit import Chem

def create_logger(name: str, task_id: int) -> logging.Logger:
    """
    Creates a logger with a stream handler and two file handlers.

    The stream handler prints to the screen depending on the value of `quiet`.
    One file handler (verbose.log) saves all logs, the other (quiet.log) only saves important info.

    :param save_dir: The directory in which to save the logs.
    :return: The logger.
    """
    logging.basicConfig(
            filemode='w+',
            level=logging.INFO)
    logger = logging.getLogger(name)
    logger.propagate = False
    file_name = f'{name}_{task_id}.log'
    try:
        os.remove(file_name)
    except:
        pass
    fh = logging.FileHandler(filename=file_name)
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)

    return logger

def write_mol_to_sdf(mol, path, confIds = [0]):
    if isinstance(confIds, int):
        confIds = [confIds]
    writer = Chem.SDWriter(path)
    for confId in confIds:
        writer.write(mol, confId=confId)
    writer.close()

def load_sdf(path, removeHs=False, sanitize=False):
    return Chem.SDMolSupplier(path, removeHs=removeHs, sanitize=sanitize)

class DoneJobsRecord(object):
    """
    class to record completed jobs
    """
    def __init__(self):
        self.FF_conf = []
        self.XTB_opt_freq = []
        self.DFT_opt_freq = []
        self.COSMO = []
        self.WFT_sp = []
        self.QM_desp = []
    
    def save(self, project_dir, args):
        with open(os.path.join(project_dir, f"done_jobs_record_{args.task_id}.yml"), "w+") as fh:
            yaml.dump(vars(self), stream=fh)

    def load(self, project_dir, args):
        with open(os.path.join(project_dir,f"done_jobs_record_{args.task_id}.yml"), "r") as fh:
            content = yaml.load(stream=fh, Loader=yaml.Loader)
        for job, molids in content.items():
            setattr(self, job, molids)

done_jobs_record = DoneJobsRecord()