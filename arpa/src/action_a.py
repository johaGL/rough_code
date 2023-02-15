import argparse

def action_a_args():
    parser = argparse.ArgumentParser(prog="arpa.src.action_a" )
                                   # , description="action_a options (A)"
    parser.add_argument('wdir', type=str,
                        help="working directory, absolute path")
    parser.add_argument('-td' ,'--todo', help='what to do', required=True)

    return parser

if __name__ == "__main__":
    parser = action_a_args()
    args = parser.parse_args()

    print('this is action a', args.todo)

