import argparse
import sys
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

parser = argparse.ArgumentParser(
                                 prog='Exercise 6', 
                                 description='Disulphide bonds search'
                                 )

# Use add_argument for all needed arguments
# https://docs.python.org/3/library/argparse.html

# Format for a optional text argument, default values can be indicated.
parser.add_argument(
    '-d', '--distance',
    type=float,
    default=2.0,
    help='Distance for a dissulfide bond, default - 2.0'
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


for at in st.get_atoms():
    if at.id == 'SG' and at.get_parent().get_resname() == 'CYS':
        select.append(at)
        print(f"ATOM: {at.get_parent().get_resname()}, {at.get_parent().id[1]}, {at.id}")

# Preparing search
if len(select) == 0:
    print("No cysteine SG atoms found in this structure — no disulphide bonds possible.")
    sys.exit(0)

nbsearch = NeighborSearch(select)

print("NBSEARCH:")

contacts = []

for a, b in nbsearch.search_all(MAXDIST):
    res_a = a.get_parent()
    res_b = b.get_parent()

    if res_a == res_b:
        continue  # same residue, not a peptide bond
    pair = (a, b) if res_a.id[1] <= res_b.id[1] else (b, a)
    contacts.append(pair)

contacts.sort(key=lambda pair: (pair[0].get_parent().id[1], pair[1].get_parent().id[1]))

ncontact = 1

for at1, at2 in contacts:
    res_1 = at1.get_parent()
    res_2 = at2.get_parent()
    print(f"Disulphide bond {ncontact}: "
          f"{res_1.get_resname()} {res_1.id[1]}  <-->  "
          f"{res_2.get_resname()} {res_2.id[1]}")
    ncontact += 1
