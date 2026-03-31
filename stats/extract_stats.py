import sys
from pathlib import Path

import os, sys; sys.path.insert(0, os.getcwd()) # allow imports from root folder
from hire2fa.pdb_table import PDBTable
from hire2fa.residue import Residue
from hire2fa.geometry import Geometry
from hire2fa.mapping import \
    BASENAMES_FOR_CENTROID_1, \
    BASENAMES_FOR_CENTROID_2, \
    NAME_CENTROID as NC


# ------------------------------------------------------------------------------
def calc_base_centroid(residue: Residue, basenames: str, i: int) -> PDBTable:
    data = residue.get_functional_group(basenames[residue.restype])
    if not len(data): return data # empty table if no atoms of the base are present

    coords = data.get_coords_array()
    centroid = Geometry.calc_center_of_mass(coords, data.sections["atomname"])
    out = {k:(tup[0],) for k,tup in data.sections.items()}
    out["atomname"] = (NC(i),)
    out["xcoords"]  = (f"{centroid[0]:.3f}",)
    out["ycoords"]  = (f"{centroid[1]:.3f}",)
    out["zcoords"]  = (f"{centroid[2]:.3f}",)
    return PDBTable(sections = out)


# ------------------------------------------------------------------------------
def extract_geometries(residue: Residue) -> tuple[str, str, Geometry]:
    def pack_geo(atomname: str, ref: tuple[str, str, str], alias: str):
        geo = residue.init_fa_geometries(*ref)
        return (alias, residue.restype, residue.calc_geometry_values(geo, atomname))

    geometries = (
        pack_geo("OP1", ("O3'",   "P", "O5'"), "op1"),
        pack_geo("OP2", ("O3'",   "P", "O5'"), "op2"),
        pack_geo("C5'", ("C1'", "C4'", "O5'"), "c5" ),
        pack_geo("O4'", ("C5'", "C4'", "C1'"), "o4" ),
        pack_geo("C3'", ("O5'", "C4'", "C1'"), "c3" ),
        pack_geo("C2'", ("C4'", "O4'", "C1'"), "c2" ),
        pack_geo("O2'", ("C4'", "C1'", "C2'"), "o2" ),
    )
    match residue.restype:
        case 'U': geometries += (
            pack_geo("C6", (  "P", "C1'", NC(1)), "buC6"),
            pack_geo("N1", ("C1'", NC(1),  "C6"), "buN1"),
            pack_geo("C2", ( "C6",  "N1", NC(1)), "buC2"),
            pack_geo("O2", ( "C6",  "N1", NC(1)), "buO2"),
            pack_geo("N3", ( "C6",  "N1", NC(1)), "buN3"),
            pack_geo("C4", ( "N1",  "C6", NC(1)), "buC4"),
            pack_geo("O4", (NC(1),  "N3",  "C4"), "buO4"),
            pack_geo("C5", ( "N1",  "C6", NC(1)), "buC5"),
        )
        case 'C': geometries += (
            pack_geo("C6", ("C4'", "C1'", NC(1)), "bcC6"),
            pack_geo("N1", ("C1'", NC(1),  "C6"), "bcN1"),
            pack_geo("C2", ( "C6",  "N1", NC(1)), "bcC2"),
            pack_geo("O2", ( "C6",  "N1", NC(1)), "bcO2"),
            pack_geo("N3", ( "C6",  "N1", NC(1)), "bcN3"),
            pack_geo("C4", ( "N1",  "C6", NC(1)), "bcC4"),
            pack_geo("N4", (NC(1),  "N3",  "C4"), "bcN4"),
            pack_geo("C5", ( "N1",  "C6", NC(1)), "bcC5"),
        )
        case 'A': geometries += (
            pack_geo("N9", ("C1'", NC(2), NC(1)), "baN9"),
            pack_geo("C4", ( "N9", NC(2), NC(1)), "baC4"),
            pack_geo("N3", ( "N9", NC(1), NC(2)), "baN3"),
            pack_geo("C2", ( "N9", NC(1), NC(2)), "baC2"),
            pack_geo("N1", ( "N9", NC(1), NC(2)), "baN1"),
            pack_geo("C6", ( "N9", NC(1), NC(2)), "baC6"),
            pack_geo("N6", ( "N9", NC(1), NC(2)), "baN6"),
            pack_geo("C5", ( "N9", NC(1), NC(2)), "baC5"),
            pack_geo("N7", ( "N9", NC(2), NC(1)), "baN7"),
            pack_geo("C8", (NC(2), NC(1),  "N9"), "baC8"),
        )
        case 'G': geometries += (
            pack_geo("N9", ("C1'", NC(2), NC(1)), "bgN9"),
            pack_geo("C4", ( "N9", NC(2), NC(1)), "bgC4"),
            pack_geo("N3", ( "N9", NC(1), NC(2)), "bgN3"),
            pack_geo("C2", ( "N9", NC(1), NC(2)), "bgC2"),
            pack_geo("N2", ( "N9", NC(1), NC(2)), "bgN2"),
            pack_geo("N1", ( "N9", NC(1), NC(2)), "bgN1"),
            pack_geo("C6", ( "N9", NC(1), NC(2)), "bgC6"),
            pack_geo("O6", ( "N9", NC(1), NC(2)), "bgO6"),
            pack_geo("C5", ( "N9", NC(1), NC(2)), "bgC5"),
            pack_geo("N7", ( "N9", NC(2), NC(1)), "bgN7"),
            pack_geo("C8", (NC(2), NC(1),  "N9"), "bgC8"),
        )

    return filter(lambda tup: tup[2] is not None, geometries)


# ------------------------------------------------------------------------------
def process_pdb(path_pdb: str | Path) -> None:
    pdb = PDBTable.read_pdb(path_pdb)
    chains = Residue.get_pdb_chains(pdb, init_with_cg_data = False)
    residues = [residue for chain in chains for residue in chain]

    for residue in residues:
        residue.fa_data.append_table(calc_base_centroid(residue, BASENAMES_FOR_CENTROID_1, 1))
        if residue.restype not in ("A", "G"): continue
        residue.fa_data.append_table(calc_base_centroid(residue, BASENAMES_FOR_CENTROID_2, 2))

    Residue.apply_o3_shift(chains, do_cg = False)

    geometries = (tup_geo for residue in residues for tup_geo in extract_geometries(residue))

    with open(PATH_CSV_STATS, 'a') as file:
        file.write('\n'.join(
            ','.join((
                path_pdb.stem, restype, atomname,
                f"{geo.dist:.3f}", f"{geo.angle:.6f}", f"{geo.dihed:.6f}"
            ))
            for atomname,restype,geo in geometries
        ) + '\n')


# ------------------------------------------------------------------------------
def main():
    with open(PATH_CSV_STATS, 'w') as file:
        file.write("pdb,resname,atomname,dist,angle,dihed\n")

    for path_pdb in FOLDER_PDB_FA.glob("*.pdb"):
        print(">>> Processing", path_pdb.name)
        process_pdb(path_pdb)


################################################################################
if __name__ == "__main__":
    FOLDER_PDB_FA  = Path(sys.argv[1])
    PATH_CSV_STATS = Path(sys.argv[2])
    PATH_CSV_STATS.parent.mkdir(parents = True, exist_ok = True)
    main()


################################################################################
# python3 stats/extract_stats.py testdata/dataset_fa testdata/stats.csv
