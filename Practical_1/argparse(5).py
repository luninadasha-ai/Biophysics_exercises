import argparse
import sys
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

parser = argparse.ArgumentParser(
                                 prog='Exercise 5', 
                                 description='Generate a list of backbone connectivity (i.e. which residues are linked by ordinary peptide bonds).'
                                 )

# Use add_argument for all needed arguments
# https://docs.python.org/3/library/argparse.html

# Format for a optional text argument, default values can be indicated.
parser.add_argument(
    '-d', '--distance',
    type=float,
    default=2.5,
    help='Distance for a peptide, default - 2.5'
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
Dist = args.distance
Pdb_file = args.PDB_file

print(Dist, Pdb_file) 

MAXDIST = Dist  # Define distance for a  contact

PDB_parser = PDBParser()

# load structure from PDB file

st = PDB_parser.get_structure(args.PDB_file, args.PDB_file)

select = []

elems = ['C', 'N']

for at in st.get_atoms():
    if at.id in elems:
        select.append(at)
        print(f"ATOM: {at.get_parent().get_resname()}, {at.get_parent().id[1]}, {at.id}")

# Preparing search
nbsearch = NeighborSearch(select)

print("NBSEARCH:")

contacts = []

for a, b in nbsearch.search_all(MAXDIST):
    if a.id == 'C' and b.id == 'N':
        c_atom, n_atom = a, b
    elif a.id == 'N' and b.id == 'C':
        c_atom, n_atom = b, a
    else:
        continue  # not a C-N pair

    res_c = c_atom.get_parent()
    res_n = n_atom.get_parent()

    if res_c == res_n:
        continue  # same residue, not a peptide bond

    if res_c.get_parent().id != res_n.get_parent().id:
        continue  # different chains, not a peptide bond

    contacts.append((c_atom, n_atom))

contacts.sort(key=lambda pair: pair[0].get_parent().id[1])

ncontact = 1

for c_atom, n_atom in contacts:
    res_c = c_atom.get_parent()
    res_n = n_atom.get_parent()
    print(f"Bond {ncontact}: "
          f"{res_c.get_resname()} {res_c.id[1]} (C)  <-->  "
          f"{res_n.get_resname()} {res_n.id[1]} (N)")
    ncontact += 1
