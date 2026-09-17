import argparse
import sys
from Bio.PDB.PDBParser import PDBParser

parser = argparse.ArgumentParser(
                                 prog='Exercise 2', 
                                 description='Generate a list of all atoms for a given residue number'
                                 )

# Use add_argument for all needed arguments
# https://docs.python.org/3/library/argparse.html

# Format for a optional text argument, default values can be indicated.
parser.add_argument('residue_num',
                    help='Residue number, optionally with chain as CHAIN:NUM (e.g. A:35)'
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
Residue = args.residue_num
Pdb_file = args.PDB_file

print(Residue, Pdb_file)

if ':' in Residue:
    chain_id, resnum_str = Residue.split(':')
else:
    chain_id = None
    resnum_str = Residue

try:
    resnum = int(resnum_str)
except ValueError:
    sys.exit(f"Error: '{resnum_str}' is not a valid residue number")

PDB_parser = PDBParser()

st = PDB_parser.get_structure(args.PDB_file, args.PDB_file)

selected = []

for at in st.get_atoms():
    residue = at.get_parent()
    chain = residue.get_parent()
    if residue.id[1] == resnum and (chain_id is None or chain.id == chain_id):
        selected.append(at)

selected.sort(key=lambda atom: atom.get_serial_number())

print("Coordinates:")
for atom in selected:
    residue = atom.get_parent()
    print(f"{residue.get_resname()}, {residue.id}, {atom.get_name()}, {atom.get_coord()}")
