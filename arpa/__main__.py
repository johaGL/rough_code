"""
in a README must state:
HELP:
  to see all options across each module at once, run
$ python -m arpa
without options, or alternative, launch help independently:
$ python -m arpa.src.action_a --help


"""
from .src.action_b import action_b_Arguments
from .src.action_a import action_a_args

if __name__ == "__main__":

    print("\n * Options across all DIMet modules * \n")

    parserB = action_b_Arguments()
    parserB.print_help()
    print("")
    parserA = action_a_args()
    parserA.print_help()




