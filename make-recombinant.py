from Bio import SeqIO
import json
import argparse
import sys
import re

parser = argparse.ArgumentParser("Generate recombinant from examples")
parser.add_argument("terms", nargs='+',
                    help="Specify a recombinant with a space-delimited sequence "
                         "of terms, e.g., python3 make-recombiant.py BA.1 10000 BA.2 20000 BA.3")
args = parser.parse_args()


def get_coord(mutation):
    position = re.sub("[^0-9]+([0-9]+)[^0-9]+", "\\1", mutation['mutation'])
    if position.isdigit():
        return int(position)
    return None


# load virus definitions
vprop = json.load(open("virus_properties.json"))

"""
>>> vprop.keys()
dict_keys(['schemaVersion', 'comment', 'variants'])
>>> vprop['variants'][0].keys()
dict_keys(['NextstrainClade', 'PangoLineage', 'Letter', 'WhoLabel', 'Other',
           'WhoClass', 'mutations', 'name'])
>>> vprop['variants'][0]['mutations'][0]
{'mutation': 'C23709T', 'proportion': 0.9836499757437491, 'count': 525154}
"""


variants = {}
all_coords = set()
for variant in vprop['variants']:
    variants.update({variant['PangoLineage']: variant})
    for mutation in variant['mutations']:
        coord = get_coord(mutation)
        all_coords.update({coord})

all_coords = list(all_coords)
all_coords.sort()


# check tokens
if len(args.terms) % 2 != 1:
    print("ERROR: You must enter an odd number of terms")
    sys.exit()

parents = []
breakpoints = []

for i in range(len(args.terms)):
    if i % 2 == 0:
        pango = args.terms[i]
        if pango not in variants:
            print(f"ERROR: Unrecognized PANGO lineage {pango}: please choose from:")
            print(f"{', '.join([p for p in variants.keys() if p])}")
            sys.exit()
        parents.append(pango)
    else:
        pos = args.terms[i]
        if not pos.isdigit():
            print(f"ERROR: Term {pos} must be a non-negative integer")
            sys.exit()
        pos = int(pos)
        if pos < 0 or pos > 29903:
            print(f"ERROR: Position {pos} must be between 0 and 29903")
            sys.exit()
        breakpoints.append(pos)

# load reference genome
record = SeqIO.read(open("reference.fasta"), 'fasta')
seqlist = list(record.seq)  # make mutable

for variant in parents:
    mutations = variants[variant]['mutations']
    intermed = [(get_coord(m), m['mutation']) for m in mutations]
    intermed.sort()
    for coord, m in intermed:
        if breakpoints and coord > breakpoints[0]:
            breakpoints.pop(0)
            break

        # apply mutation to sequence (coordinates are 1-index)
        seqlist[coord-1] = m[-1]

# write recombinant sequence to console
print(''.join(seqlist))
