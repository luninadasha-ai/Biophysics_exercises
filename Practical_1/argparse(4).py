import argparse
import sys
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

parser = argparse.ArgumentParser(
                                 prog='Exercise 4', 
                                 description='Generate a list of all CA atoms of given residue type with coordinates'
                                 )

# Use add_argument for all needed arguments
# https://docs.python.org/3/library/argparse.html

# Format for a optional text argument, default values can be indicated.
parser.add_argument('residue_type',
                    help='Residue type'
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
Residue_type = args.residue_type
Pdb_file = args.PDB_file

print(Residue_type, Pdb_file)

one_to_three = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}
 
user_input = args.residue_type.strip().upper()
 
if len(user_input) == 1:
    if user_input not in one_to_three:
        sys.exit(f"Error: '{user_input}' is not a valid one-letter amino acid code")
    Residue_type = one_to_three[user_input]
elif len(user_input) == 3:
    if user_input not in one_to_three.values():
        sys.exit(f"Error: '{user_input}' is not a valid three-letter amino acid code")
    Residue_type = user_input
else:
    sys.exit(f"Error: '{user_input}' is not a valid residue code (expected 1 or 3 letters)")
 
print(Residue_type, Pdb_file)

PDB_parser = PDBParser()

# load structure from PDB file

st = PDB_parser.get_structure(args.PDB_file, args.PDB_file)

select = []

#Select only CA atoms

for at in st.get_atoms():
    if at.id == 'CA' and at.get_parent().get_resname() == Residue_type :
        select.append(at)

select.sort(key=lambda atom: atom.get_serial_number())
for at in select:
    print(f"ATOM: {at.get_parent().get_resname()}, {at.coord}, {at.id}")

