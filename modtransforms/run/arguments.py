from os import path
from argparse import ArgumentParser

def parse_arguments():
    parser = ArgumentParser(prog=path.splitext(path.basename(__file__))[0],
                            description='Update positions using MOD File')
    parser.add_argument('--mod', required=True, help='Required MOD File')
    args = parser.parse_args()
