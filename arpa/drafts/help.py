"""
this file works in arpa/ location, not in arpa/drafts/
this is ugly, beter see __main__ and __init__
"""


import subprocess
import os


if __name__ == "__main__":
    os.chdir(os.path.dirname(__file__))
    print("")
    subprocess.run(["python","-m", "src.action_a", "--help"])
    print("")
    subprocess.run(["python", "-m", "src.action_b", "--help"])
