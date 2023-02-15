import argparse

def action_b_Arguments():
    parser = argparse.ArgumentParser(prog="python -m arpa.src.action_b")
    parser.add_argument('wdir', type=str,
                        help="working directory, absolute path")
    parser.add_argument('config', type=str, help="xxxxxx, sdssss")
    parser.add_argument('-e' ,'--eat', help='what to eat')

    return parser

if __name__ == "__main__":
    parser = action_b_Arguments()
    args = parser.parse_args()
    print('this is action b !', args.eat)
    print("wdir :", args.wdir, "\nconfig :", args.config)
