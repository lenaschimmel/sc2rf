import csv
from typing import NamedTuple
from termcolor import colored, cprint
import fileinput
import json
import sys

from xarray import Coordinate

colors = ['red', 'green', 'blue', 'yellow', 'magenta', 'cyan']

# I removed "ORF" from the names, because often we only see the first one or two letters of a name, and "ORF" providex no information
genes = {
    '1a': (  266, 13468),
    '1b': (13468, 21555),
    'S':  (21563, 25384),
    'E':  (26245, 26472),
    '3a': (25393, 26220),
    'M':  (26523, 27191),
    '6':  (27202, 27387),
    '7a': (27394, 27759),
    '7b': (27756, 27887),
    '8':  (27894, 28259),
    'N':  (28274, 28283), # my algorithm does not like that ORF9b is inside N, so I split N in two halves
    '9b': (28284, 28577),
    'N':  (28578, 29533),
    '': (29534, 99999), # probably nothing, but need something to mark the end of N
}

def main():
    global reference
    global mappings
    print("Reading name mappings, reference genome, lineage definitions...")
    mappings = read_mappings('mapping.csv')
    reference = read_fasta('reference.fasta')['MN908947 (Wuhan-Hu-1/2019)']
    all_examples = read_examples('virus_properties.json')

    print("Done.\nReading actual input.")s
    all_samples = read_subs_from_fasta('sample.fasta')
    print("Done.")

    # The current algorithm relies on "unique mutations", that is, those which only occur in one of
    # the example clades. Having to many similar clades reduces the number of unique mutations
    # too much to be useful.
    # As a quick fix, let's filter the examples to a much smaller list:
    used_clades = ['20I','20H','20J','21A','21K','21L']
    used_examples = dict()
    for ex_name, ex in all_examples.items():
        if ex_name in used_clades:
            used_examples[ex_name] = ex

    calculate_relations(used_examples)

    match_sets = dict()

    print("Scanning input for matches against linege definitons...")
    for sa_name, sa in all_samples.items():
        matching_example_names = []
        for ex_name, ex in used_examples.items():
            matches_count = len(sa['subs_set'] & ex['unique_subs_set'])
            #print(f"    {matches_count}")
            #matches_percent = int(matches_count / len(ex['unique_subs_set']) * 100)
            if matches_count > 0: # theoretically > 0 already gives us recombinants, but they are much more likely to be errors or coincidences
                matching_example_names.append(ex_name)# f"  Unique {ex_name}: {matches_percent}% ({matches_count} of {len(ex['unique_subs_set'])})")

        matching_examples_tup = tuple(matching_example_names)

        if len(matching_example_names) > 1:
            #print(f"{sa_name} is a possible recombinant of {len(matching_example_names)} lineages: {matching_example_names}")
            if match_sets.get(matching_examples_tup):
                match_sets[matching_examples_tup].append(sa)
            else:
                match_sets[matching_examples_tup] = [sa]

    print("Done.\nPriniting detailed analysis:\n\n")

    for example_names, samples in match_sets.items():
        show_matches(used_examples, example_names, samples)

def pretty_name(clade_name):
    global mappings
    clade = mappings['by_clade'].get(clade_name)
    if clade:
        label = clade.get('WhoLabel')
        pango = clade.get('PangoLineage')
        if label and pango:
            return f"{label} ({pango} / {clade_name})"
        elif pango:
            return f"{pango} / {clade_name}"
        else:
            return clade_name
    else:
        return f"Unknown ({clade_name})"

def read_examples(path):
    with open(path, newline='') as jsonfile:
        props = json.load(jsonfile)
        assert props['schemaVersion'] == '1.10.0'
        sequences = {}
        for clade, subs in props['nucMutLabelMapReverse'].items():
            name = clade # TODO get pango-names
            subs_dict = dict()
            for s in subs:
                s = s.strip()
                if len(s) > 0:
                    sub = parse_short_sub(s)
                    subs_dict[sub.coordinate] = sub

            #print(f"Read {name} with subs: {subs_dict}")

            sequences[name] = {
                'clade': clade,
                'name': name,
                'subs_dict': subs_dict,
                'subs_list': list(subs_dict.values()),
                'subs_set': set(subs_dict.values()),
                'missings': []
            }
        return sequences

def read_fasta(path):
    sequences = dict()
    with open(path, newline='') as fasta:
        for line in fasta:
            if(line[0] == '>'):
                current_name = line[1:].strip()
                sequences[current_name] = ""
            else:
                sequences[current_name] += line.strip()

    return sequences;

def read_subs_from_fasta(path):
    fastas = read_fasta(path)
    sequences = dict()
    start_n = -1
    for name, fasta in fastas.items():
        subs_dict = dict()
        missings = list()
        if len(fasta) != len(reference):
            print(f"Sequence {name} not properly aligned, length is {len(fasta)} instead of {len(reference)}.")
        else:
            for i in range(1, len(reference) + 1):
                r = reference[i - 1]
                s = fasta[i - 1]
                if s == 'N':
                    if start_n == -1:
                        start_n = i
                elif start_n >= 0:
                    missings.append((start_n, i - 1))
                    start_n = -1
                    
                if s != 'N' and r != s:
                    subs_dict[i] = Sub(r, i, s)
                    

            sequences[name] = {
                'name': name,
                'subs_dict': subs_dict,
                'subs_list': list(subs_dict.values()),
                'subs_set': set(subs_dict.values()),
                'missings': missings
            }

    return sequences

class Sub(NamedTuple):
    ref: str
    coordinate: int
    mut: str

def parse_sub(s):
    return Sub(s[0], int(s[1:-1]), s[-1])

def parse_short_sub(s):
    coordinate = int(s[0:-1])
    return Sub(reference[coordinate], coordinate, s[-1])


def prunt(s, color=None):
    if color:
        cprint(s, color, end="")
    else:
        print(s, end="")

def fixed_len(s, l):
    trunc = s[0:l]
    return trunc.ljust(l)

def show_matches(all_examples, example_names, samples):
    ml = 52

    examples = [all_examples[name] for name in example_names]

    coords = set()
    refs = dict() 
    for ex in examples:
        for sub in ex['subs_list']:
            coords.add(sub.coordinate)
            refs[sub.coordinate] = sub.ref;
   
    # show all sample mutations:
    # for sa in samples:
    #     for sub in sa['subs_list']:
    #         coords.add(sub.coordinate)
    #         refs[sub.coordinate] = sub.ref;
    
    ordered_coords = list(coords)
    ordered_coords.sort()

    print(f"\n\nPotential recombinants between {example_names}:\n")
   
    ###### SHOW COORDS

    for exp in range(5,0,-1):
        div = 10**(exp-1)

        if exp == 5:
            prunt(fixed_len("coordinates", ml + 1))
        else:
            prunt(' ' * (ml+1))

        for coord in ordered_coords:
            if coord//div > 0:
                prunt((coord//div)%10)
            else:
                prunt(' ')
            #print(f"{coord} // {div} = {(coord//div)}")
        print()
    print()

    ###### SHOW GENES
    prunt(fixed_len("genes", ml + 1))

    current_name = ''
    color_index = 0
    current_color = colors[color_index]
    text_index = 0

    for coord in ordered_coords:
        for name, limits in genes.items():
            if coord >= limits[0] and coord <= limits[1]:
                if current_name != name:
                    current_name = name
                    color_index = (color_index + 1) % len(colors)
                    current_color = colors[color_index]
                    text_index = 0
        char = ' '
        if len(current_name) > text_index:
            char = current_name[text_index]
        cprint(char, current_color, None, attrs=['reverse'], end='')
        text_index += 1
    print()

    ###### SHOW REF
    
    prunt(fixed_len("ref", ml + 1))
    for coord in ordered_coords:
        prunt(refs[coord])
    print()
    print()

    ###### SHOW EXAMPLES

    color_index = 0
    color_by_name = dict()

    for ex in examples:
        current_color = colors[color_index]
        color_by_name[ex['name']] = current_color
        prunt(fixed_len(pretty_name(ex['name']), ml) + ' ', current_color)
        for coord in ordered_coords:
            if(ex['subs_dict'].get(coord)):
                prunt(ex['subs_dict'][coord].mut, current_color)
            else:
                prunt(".")
        print()
        color_index += 1
    print()

    ###### SHOW SAMPLES
    current_color = 'grey'

    for sa in samples:
        #current_color = colors[color_index]
        #color_by_name[sa['name']] = current_color

        prev_definitive_match = None
        breakpoints = 0
        definitives_since_breakpoint = 0
        definitives_count = []

        prunt(fixed_len(sa['name'], ml) + ' ', current_color)
        for coord in ordered_coords:
            if is_missing(coord, sa['missings']):
                cprint('N', 'white', attrs=['reverse'], end='')
            else:
                if(sa['subs_dict'].get(coord)): # sample has sub here
                    matching_exs = []
                    for ex in examples:
                        if ex['subs_dict'].get(coord) and ex['subs_dict'].get(coord).mut == sa['subs_dict'][coord].mut:
                            matching_exs.append(ex['name'])

                    text = sa['subs_dict'][coord].mut
                    fg = 'white'
                    bg = None
                    attrs=['bold']

                    if(len(matching_exs) == 0): # none of the examples match - private mutation
                        bg = 'on_cyan'
                    elif(len(matching_exs) == 1): # exactly one of the examples match - definite match
                        fg = color_by_name[matching_exs[0]]
                        if matching_exs[0] != prev_definitive_match:
                            if prev_definitive_match:
                                breakpoints += 1
                            if definitives_since_breakpoint:
                                definitives_count.append((prev_definitive_match,definitives_since_breakpoint))
                            prev_definitive_match = matching_exs[0]
                            definitives_since_breakpoint = 0
                        definitives_since_breakpoint += 1
                    elif(len(matching_exs) < len(examples)): # more than one, but not all examples match - can't provide proper color
                         #bg = 'on_blue'
                         attrs = ['bold', 'underline']
                    # else: all examples match

                    cprint(text, fg, bg, attrs=attrs, end='')
                else: # sample does not have sub here
                    matching_exs = []
                    for ex in examples:
                        if not ex['subs_dict'].get(coord):
                            matching_exs.append(ex['name'])

                    text = '•'
                    #text = refs[coord]
                    fg = 'white'
                    bg = None
                    attrs = []

                    if(len(matching_exs) == 0):  # none of the examples match - private reverse mutation
                        bg = 'on_magenta'
                    elif(len(matching_exs) == 1): # exactly one of the examples match - definite match
                        fg = color_by_name[matching_exs[0]]
                        if matching_exs[0] != prev_definitive_match:
                            if prev_definitive_match:
                                breakpoints += 1
                            if definitives_since_breakpoint:
                                definitives_count.append((prev_definitive_match,definitives_since_breakpoint))
                            prev_definitive_match = matching_exs[0]
                            definitives_since_breakpoint = 0
                        definitives_since_breakpoint += 1
                    elif(len(matching_exs) < len(examples)): # more than one, but not all examples match - can't provide proper color
                         #bg = 'on_yellow'
                         attrs = ['underline']
                    # else: all examples match (which should not happen, because some example must have a mutation here)
                    
                    cprint(text, fg, bg, attrs=attrs, end='')
        if definitives_since_breakpoint:
            definitives_count.append((prev_definitive_match,definitives_since_breakpoint))

        # now transform definitive streaks: every sequence like ..., X, S, Y, ... where S is a small numer into ..., (X+Y), ...

        min_streak_length = 3

        reduced = list(filter(lambda ex_count: ex_count[1] >= min_streak_length, definitives_count))
        removed = len(definitives_count) - len(reduced)
        further_reduced = []

        last_ex = reduced[0][0];
        last_count = 0
        for (ex, count) in reduced:
            if ex != last_ex:
                further_reduced.append(last_count)
                last_count = count
                last_ex = ex
            else:
                last_count += count

        if last_count:
            further_reduced.append(last_count)
        
        prunt(f"     {len(further_reduced) - 1} breakpoint(s)")
        if removed:
            prunt(f", ignored {removed} streaks < {min_streak_length}")

        if len(further_reduced) > 5 or len(further_reduced) == 1:
            print("\x0D\x1B[2K", end="") # remove the whole line which we just printed
        else:
            print()
    
    print()

def read_mappings(path):
    with open(path, newline='') as csvfile:
        mappings = { 
            'by_clade': dict()
        }
        reader = csv.DictReader(csvfile)
        line_count = 0
        for row in reader:
            mappings['by_clade'][row['NextstrainClade']] = row
    return mappings

def read_subs(path, delimiter = ',', max_lines = -1):
    with open(path, newline='') as csvfile:
        sequences = {}
        reader = csv.DictReader(csvfile, delimiter=delimiter)
        line_count = 0
        for row in reader:
            subs_dict = dict()
            missings = list()
            for s in row['substitutions'].split(","):
                s = s.strip()
                if len(s) > 0:
                    sub = parse_sub(s)
                    subs_dict[sub.coordinate] = sub

            for m in row['missing'].split(","):
                m = m.strip()
                if len(m) > 0:
                    parts = m.split('-')
                    if len(parts) == 1:
                        missings.append((int(parts[0]),int(parts[0])))
                    else:
                        missings.append((int(parts[0]),int(parts[1])))

            sequences[row['seqName']] = {
                'name': row['seqName'],
                'subs_dict': subs_dict,
                'subs_list': list(subs_dict.values()),
                'subs_set': set(subs_dict.values()),
                'missings': missings
            }
           
            line_count += 1
            if max_lines != -1 and line_count == max_lines:
                break
    return sequences

def is_missing(coordinate, missings):
    for missing in missings:
        if coordinate >= missing[0] and coordinate <= missing[1]:
            return True
    return False

def calculate_relations(examples):
    for example in examples.values():
        union = set()
        for other in examples.values():
            if other is not example:
                union = union | (other['subs_set'])
        example['unique_subs_set'] = example['subs_set'] - union
        print(f"Clade  {pretty_name(example['name'])} has {len(example['subs_set'])} mutations, of which {len(example['unique_subs_set'])} are unique.")

if __name__ == '__main__':
    main()