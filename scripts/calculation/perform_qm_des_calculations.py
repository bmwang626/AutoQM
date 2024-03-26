import logging
from argparse import ArgumentParser
from radical_workflow.calculation.utils import add_shared_arguments, add_qm_des_arguments

logging.basicConfig(level=logging.INFO)

def main(args):
    logging.info("Performing QM descriptor calculations")
    

if __name__ == "__main__":
    parser = ArgumentParser()
    parser = add_shared_arguments(parser)
    parser = add_qm_des_arguments(parser)
    args = parser.parse_args()
    main(args)