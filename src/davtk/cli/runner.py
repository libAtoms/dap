import sys

from pathlib import Path
from subprocess import call

def main():
    call([str(Path(__file__).parent.parent / "scripts" / Path(__file__).stem) ] + sys.argv[1:])
