import sys
from pathlib import Path

import hire2fa as h2f

# ------------------------------------------------------------------------------
def main():
    rec = h2f.Reconstructor(
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
# hire2fa testdata/examples/1akx.cg.pdb testdata/examples/1akx.hire2fa.pdb
# python3 -m hire2fa testdata/examples/1akx.cg.pdb testdata/examples/1akx.hire2fa.pdb
