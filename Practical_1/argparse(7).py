import argparse
import sys
from Bio.PDB.PDBParser import PDBParser
import numpy as np

parser = argparse.ArgumentParser(
                                 prog='Exercise 7', 
                                 description='Print distances between all atom pairs of two given residues'
                                 )

# Use add_argument for all needed arguments
# https://docs.python.org/3/library/argparse.html

# Format for a optional text argument, default values can be indicated.
parser.add_argument('residue_1',
                    help='Name of the first residue'
                    )

parser.add_argument('residue_2',
                    help='Name of the second residue'
                    )

# Format for a required parameter, call will fail if empty. This only stores the file name, but specifiying type=file, the file is open and contents available. 
parser.add_argument('PDB_file',
                    help='Required PDB file for the program')

# Read command line into args

args = parser.parse_args()
       
# Print the parameters that has been read 
    
print ("\nSettings\n--------")

for k, v in vars(args).items():
    print ('{:10}:'.format(k), v)

print ("\nSettings, again\n---------------")    

#print the variables once assigned
Res_1 = args.residue_1
Res_2 = args.residue_2
Pdb_file = args.PDB_file

print(Res_1, Res_2, Pdb_file) 

PDB_parser = PDBParser()

st = PDB_parser.get_structure(args.PDB_file, args.PDB_file)

def get_residue(structure, res_spec):
    """Parse 'CHAIN:NUM' or just 'NUM' and return the Biopython residue."""
    if ':' in res_spec:
        chain_id, resnum_str = res_spec.split(':')
    else:
        chain_id = next(iter(structure[0].child_dict))
        resnum_str = res_spec

    resnum = int(resnum_str)
    if chain_id not in structure[0]:
        sys.exit(f"Chain '{chain_id}' not found")

    chain = structure[0][chain_id]
    if resnum not in [res.id[1] for res in chain]:
        sys.exit(f"Residue {resnum} not found in chain {chain_id}")
    return chain[resnum]

res1 = get_residue(st, args.residue_1)
res2 = get_residue(st, args.residue_2)

print("Residue 1 is", res1.get_resname())
print("Residue 2 is", res2.get_resname())

print("\nAtom1 Atom2 dist1 dist2\n-------------------------")
pairs = []
for at1 in res1.get_atoms():      # Replace get_atoms with get_atom if you get an Error!
    for at2 in res2.get_atoms():
        dist = at2 - at1     # Direct procedure with (-) to compute distances
        vector = at2.coord - at1.coord  # Or using numpy coordinates
        distance = np.sqrt(np.sum(vector ** 2))
        pairs.append((at1, at2, dist, distance))

pairs.sort(key=lambda p: (p[0].get_serial_number(), p[1].get_serial_number()))

for at1, at2, dist, distance in pairs:
    print(at1, at2, dist, distance)

