import sys
from pathlib import Path

from reconstructor import Reconstructor

# ------------------------------------------------------------------------------
def main():
    rec = Reconstructor(
        path_pdb_cg = PATH_PDB_CG,
        path_model = Path(__file__).parent / "_data/model.json"
    )
    rec.reconstruct()
    rec.export_reconstructed(PATH_PDB_FA)


################################################################################
if __name__ == "__main__":
    PATH_PDB_CG = Path(sys.argv[1])
    PATH_PDB_FA = Path(sys.argv[2])
    main()


################################################################################
# python3 src/hire2fa/cg2fa.py testdata/conf_chainid.pdb testdata/1akx.fa.pdb
