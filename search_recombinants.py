#!/usr/bin/env python3

import csv
from typing import NamedTuple
from numpy import int0
from termcolor import colored, cprint
import fileinput
import json
import argparse

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

class Interval:
    def __init__(self, string):
        # TODO allow multiple separators, see https://stackoverflow.com/questions/1059559/split-strings-into-words-with-multiple-word-boundary-delimiters
        parts = string.split('-')
        if len(parts) == 1:
            self.min = int(parts[0])
            self.max = int(parts[0])
        elif len(parts) == 2:
            self.min = int(parts[0]) if parts[0] else None
            self.max = int(parts[1]) if parts[1] else None
        else:
            raise ValueError('invalid interval: ' + string)

    def matches(self, num):
        return num >= self.min and num <= self.max


def main():
    global mappings
    mappings = read_mappings('mapping.csv')
    clade_names = list(mappings['by_clade'].keys())

    parser = argparse.ArgumentParser(
        description='Analyse SARS-CoV-2 sequences for potential, unknown recombinant variants.', 
        epilog='An Interval can be a single number ("3"), a closed interval ("2-5" ) or an open one ("4-" or "-7").'
        ' The limts are inclusive. Only positive numbers are supported.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('input', nargs='+', help='input sequences to test, as aligned .fasta file(s)')
    parser.add_argument('--parents', '-p', default='2-4', metavar='INTERVAL', type=Interval, help='Allowed umber of potential parents of a recombinant.')
    parser.add_argument('--breakpoints', '-b', default='1-4', metavar='INTERVAL', type=Interval, help='Allowed number of breakpoints in a recombinant.')
    parser.add_argument('--clades', '-c', nargs='*', default=['20I','20H','20J','21A','21K','21L', '21BA3'], choices=(['all'] + clade_names), help='List of clades which are considered as potential parents. Use Nextclade names, i.e. "21A". Also accepts "all".')
    parser.add_argument('--unique', '-u', default=2, type=int,  metavar='NUM', help='Minimum of substitutions in a sample which are unique to a potential parent clade, so that the clade will be considered.')
    parser.add_argument('--max-intermission-length', '-l',  metavar='NUM', default=2, type=int, help='The maximum length of an intermission in consecutive substitutions. Intermissions are stretches to be ignored when counting breakpoints.')
    parser.add_argument('--max-intermission-count', '-i',  metavar='NUM', default=8, type=int, help='The maximum number of intermissions which will be ignored. Surplus intermissions count towards the number of breakpoints.')
    parser.add_argument('--max-name-length', '-n',  metavar='NUM', default=30, type=int, help='Only show up to NUM characters of sample names.')
    parser.add_argument('--max-ambiguous', '-a',  metavar='NUM', default=50, type=int, help='Maximum number of ambiguous nucs in a sample before it gets ignored.')
    parser.add_argument('--force-all-parents', '-f', action='store_true', help='Force to consider all clades as potential parents for all sequences. Only useful for debugging.')
    parser.add_argument('--select-sequences', '-s', default='0-999999', metavar='INTERVAL', type=Interval, help='Use only a specific range of inpur sequences. DOES NOT YET WORK WITH MULTIPLE INPUT FILES.')

    global args
    args = parser.parse_args()

    used_clades = args.clades
    if used_clades == ['all']:
        used_clades = clade_names

    if args.force_all_parents and not args.parents.match(len(used_clades)):
        print("The number of allowed parents, the number of selected clades and the --force-all-parents conflict so that the results must be empty.")
        return

    global reference
    print("Reading reference genome, lineage definitions...")
    reference = read_fasta('reference.fasta')['MN908947 (Wuhan-Hu-1/2019)']
    all_examples = read_examples('virus_properties.json')

    print("Done.\nReading actual input.")
    all_samples = dict()
    for path in args.input:
        read_samples = read_subs_from_fasta(path)
        all_samples = all_samples | read_samples
    print("Done.")



    used_examples = dict()
    for ex_name, ex in all_examples.items():
        if ex_name in used_clades:
            used_examples[ex_name] = ex

    calculate_relations(used_examples)

    match_sets = dict()

    print("Scanning input for matches against linege definitons...")
    for sa_name, sa in all_samples.items():
        matching_example_names = []
        if args.force_all_parents:
            matching_example_names = used_examples.keys()
        else:    
            for ex_name, ex in used_examples.items():
                matches_count = len(sa['subs_set'] & ex['unique_subs_set'])
                matches_percent = int(matches_count / len(ex['unique_subs_set']) * 100)
                #print(f"  Unique {ex_name}: {matches_percent}% ({matches_count} of {len(ex['unique_subs_set'])})")
                if matches_count >= args.unique: # theoretically > 0 already gives us recombinants, but they are much more likely to be errors or coincidences
                    matching_example_names.append(ex_name)# 

        matching_examples_tup = tuple(matching_example_names)

        if args.parents.matches(len(matching_example_names)):
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
    removed_due_to_ambig = 0
    index = 1
    for name, fasta in fastas.items():
        if args.select_sequences.matches(index):
            subs_dict = dict()
            missings = list()
            if len(fasta) != len(reference):
                print(f"Sequence {name} not properly aligned, length is {len(fasta)} instead of {len(reference)}.")
            else:
                ambiguous_count = 0
                for i in range(1, len(reference) + 1):
                    r = reference[i - 1]
                    s = fasta[i - 1]
                    if s == 'N' or s == '-':
                        if start_n == -1:
                            start_n = i
                    elif start_n >= 0:
                        missings.append((start_n, i - 1))
                        start_n = -1
                        
                    if s != 'N' and s != '-' and r != s:
                        subs_dict[i] = Sub(r, i, s)

                    if not s in "AGTCN-":
                        ambiguous_count += 1

                if ambiguous_count <= args.max_ambiguous:
                    sequences[name] = {
                        'name': name,
                        'subs_dict': subs_dict,
                        'subs_list': list(subs_dict.values()),
                        'subs_set': set(subs_dict.values()),
                        'missings': missings
                    }
                else:
                    removed_due_to_ambig += 1
        index += 1
    if removed_due_to_ambig:
        print(f"Removed {removed_due_to_ambig} of {len(fastas)} sequences with more than { args.max_ambiguous} ambiguous nucs.")

    return sequences

class Sub(NamedTuple):
    ref: str
    coordinate: int
    mut: str

def parse_sub(s):
    return Sub(s[0], int(s[1:-1]), s[-1])

def parse_short_sub(s):
    coordinate = int(s[0:-1])
    return Sub(reference[coordinate-1], coordinate, s[-1])


def prunt(s, color=None):
    if color:
        cprint(s, color, end="")
    else:
        print(s, end="")

def fixed_len(s, l):
    trunc = s[0:l]
    return trunc.ljust(l)

def show_matches(all_examples, example_names, samples):
    ml = args.max_name_length

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

    color_by_name = dict()
    color_index = 0
    for ex in examples:
        color_by_name[ex['name']] = get_color(color_index)
        color_index += 1

    ###### SHOW SAMPLES
    current_color = 'grey'
    collected_outputs = []

    for sa in samples:
        #current_color = get_color(color_index)
        #color_by_name[sa['name']] = current_color

        prev_definitive_match = None
        breakpoints = 0
        definitives_since_breakpoint = 0
        definitives_count = []

        output = ''

        output += colored(fixed_len(sa['name'], ml) + ' ', current_color)
        for coord in ordered_coords:
            if is_missing(coord, sa['missings']):
                output += colored('N', 'white', attrs=['reverse'])
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

                    output += colored(text, fg, bg, attrs=attrs)
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
                    
                    output += colored(text, fg, bg, attrs=attrs)
        if definitives_since_breakpoint:
            definitives_count.append((prev_definitive_match,definitives_since_breakpoint))

        # now transform definitive streaks: every sequence like ..., X, S, Y, ... where S is a small numer into ..., (X+Y), ...

        reduced = list(filter(lambda ex_count: ex_count[1] > args.max_intermission_length, definitives_count))
        num_intermissions = len(definitives_count) - len(reduced)
        further_reduced = []

        if len(reduced):
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

        postfix = ''
        num_breakpoints = len(further_reduced) - 1
        if num_intermissions > args.max_intermission_count:
            postfix = ' of ' + str(num_intermissions)
            num_breakpoints += (num_intermissions - args.max_intermission_count) * 2
            num_intermissions = args.max_intermission_count

        output += f"     {num_breakpoints} breakpoint(s)"
        if num_intermissions:
            output += f", ignored {num_intermissions}{postfix} intermissions <= {args.max_intermission_length}"

        if args.breakpoints.matches(num_breakpoints):
            collected_outputs.append(output)
    
    if len(collected_outputs):

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
        current_color = get_color(color_index)
        text_index = 0

        for coord in ordered_coords:
            for name, limits in genes.items():
                if coord >= limits[0] and coord <= limits[1]:
                    if current_name != name:
                        current_name = name
                        color_index += 1
                        current_color = get_color(color_index)
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

        for ex in examples:
            current_color = color_by_name[ex['name']]
            prunt(fixed_len(pretty_name(ex['name']), ml) + ' ', current_color)
            for coord in ordered_coords:
                if(ex['subs_dict'].get(coord)):
                    prunt(ex['subs_dict'][coord].mut, current_color)
                else:
                    prunt(".")
            print()
        print()

        for output in collected_outputs:
            print(output)

    print()

def get_color(color_index): 
    return colors[color_index % len(colors)]

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