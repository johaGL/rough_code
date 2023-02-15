import argparse


parser = argparse.ArgumentParser(prog="python -m src.c")
parser.add_argument('config', type=str, help="xxxxxx, sdssss")
parser.add_argument('-e' ,'--eat', help='what to eat')
args = parser.parse_args()


