import argparse
import sys
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

parser = argparse.ArgumentParser(
                                 prog='Exercise 3', 
                                 description='Determine all possible hydrogen bonds (Polar atoms at less than 3.5 Å)'
                                 )

# Use add_argument for all needed arguments
# https://docs.python.org/3/library/argparse.html

# Format for a optional text argument, default values can be indicated.
parser.add_argument(
    '-d', '--distance',
    type=float,
    default=3.5,
    help='Distance for a hydrogen bond, default - 3.5'
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

polar_elems = ['O', 'N', 'S']

for at in st.get_atoms():
    if at.element in polar_elems: 
        select.append(at)
        print(f"ATOM: {at.get_parent().get_resname()}, {at.get_parent().id[1]}, {at.id}")

# Preparing search
nbsearch = NeighborSearch(select)

print("NBSEARCH:")

contacts = [
    (a, b) if a.get_parent().id[1] <= b.get_parent().id[1] else (b, a)
    for a, b in nbsearch.search_all(MAXDIST)
    if a.get_parent() != b.get_parent()
]
contacts.sort(key=lambda pair: (pair[0].get_parent().id[1], pair[1].get_parent().id[1]))


ncontact = 1

for at1, at2 in contacts:
    print(f"Contact: {ncontact}")
    print(f"at1: {at1}, {at1.get_serial_number()}, {at1.get_parent().get_resname()}")
    print(f"at2: {at2}, {at2.get_serial_number()}, {at2.get_parent().get_resname()}")
    print()
    ncontact += 1
