#!/usr/bin/env python3

import csv
from typing import NamedTuple
from termcolor import colored, cprint
import json
import argparse
import os
import requests

colors = ['red', 'green', 'blue', 'yellow', 'magenta', 'cyan']

width_override = None

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
        self.original_string = string
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
        if self.min and num < self.min:
            return False
        if self.max and num > self.max:
            return False
            
        return True

    def __str__(self):
        return self.original_string


def main():
    global mappings
    global width_override

    # This strange line should enable handling of 
    # ANSI / VT 100 codes in windows terminal
    # See https://stackoverflow.com/a/39675059/39946
    os.system('')

    mappings = read_mappings('mapping.csv')
    clade_names = list(mappings['by_clade'].keys())

    parser = argparse.ArgumentParser(
        description='Analyse SARS-CoV-2 sequences for potential, unknown recombinant variants.', 
        epilog='An Interval can be a single number ("3"), a closed interval ("2-5" ) or an open one ("4-" or "-7").'
        ' The limits are inclusive. Only positive numbers are supported.',
        formatter_class=ArgumentAdvancedDefaultsHelpFormatter
    )
    parser.add_argument('input', nargs='*', help='input sequence(s) to test, as aligned .fasta file(s)')
    parser.add_argument('--primers', nargs='*',  metavar='PRIMER', help='Filenames of primer set(s) to visualize. The .bed formats for ARTIC and EasySeq are recognized and supported.')
    parser.add_argument('--primer-intervals', nargs='*',  metavar='INTERVAL', type=Interval, help='Coordinate intervals in which to visualize primers.')
    parser.add_argument('--parents', '-p', default='2-4', metavar='INTERVAL', type=Interval, help='Allowed number of potential parents of a recombinant.')
    parser.add_argument('--breakpoints', '-b', default='1-4', metavar='INTERVAL', type=Interval, help='Allowed number of breakpoints in a recombinant.')
    parser.add_argument('--clades', '-c', nargs='*', default=['20I','20H','20J','21A', '21K','21L', '21BA3'], help='List of clades which are considered as potential parents. Use Nextclade names, i.e. "21A". Also accepts "all".')
    parser.add_argument('--unique', '-u', default=2, type=int,  metavar='NUM', help='Minimum of substitutions in a sample which are unique to a potential parent clade, so that the clade will be considered.')
    parser.add_argument('--max-intermission-length', '-l', metavar='NUM', default=2, type=int, help='The maximum length of an intermission in consecutive substitutions. Intermissions are stretches to be ignored when counting breakpoints.')
    parser.add_argument('--max-intermission-count', '-i', metavar='NUM', default=8, type=int, help='The maximum number of intermissions which will be ignored. Surplus intermissions count towards the number of breakpoints.')
    parser.add_argument('--max-name-length', '-n', metavar='NUM', default=30, type=int, help='Only show up to NUM characters of sample names.')
    parser.add_argument('--max-ambiguous', '-a', metavar='NUM', default=50, type=int, help='Maximum number of ambiguous nucs in a sample before it gets ignored.')
    parser.add_argument('--force-all-parents', '-f', action='store_true', help='Force to consider all clades as potential parents for all sequences. Only useful for debugging.')
    parser.add_argument('--select-sequences', '-s', default='0-999999', metavar='INTERVAL', type=Interval, help='Use only a specific range of input sequences. DOES NOT YET WORK WITH MULTIPLE INPUT FILES.')
    parser.add_argument('--enable-deletions', '-d', action='store_true', help='Include deletions in lineage comparision.')
    parser.add_argument('--show-private-mutations', action='store_true', help='Display mutations which are not in any of the potential parental clades.')
    parser.add_argument('--rebuild-examples', '-r', action='store_true', help='Rebuild the mutations in examples by querying cov-spectrum.org.')
    parser.add_argument('--mutation-threshold', '-t', metavar='NUM', default=0.75, type=float, help='Consider mutations with a prevalence of at least NUM as mandatory for a clade (range 0.05 - 1.0, default: %(default)s).')
    parser.add_argument('--add-spaces', metavar='NUM', nargs='?', default=0, const=5, type=int, help='Add spaces between every N colums, which makes it easier to keep your eye at a fixed place.')
    parser.add_argument('--sort-by-id', metavar='NUM', nargs='?', default=0, const=999, type=int, help='Sort the input sequences by the ID. If you provide NUM, only the first NUM characters are considered. Useful if this correlates with meaning full meta information, e.g. the sequencing lab.')
    parser.add_argument('--verbose', '-v', action='store_true', help='Print some more information, mostly useful for debugging.')
    parser.add_argument('--update-readme', action='store_true', help=argparse.SUPPRESS)
    
    global args
    args = parser.parse_args()

    if args.update_readme:
        update_readme(parser)
        print("Readme was updated. Program exits.")
        return

    if args.rebuild_examples:
        rebuild_examples()
        if len(args.input) == 0:
            print("Examples were rebuilt, and no input sequences were provided. Program exits.")
            return
    elif len(args.input) == 0:
        print("Input sequences must be provided, except when rebuilding the examples. Use --help for more info. Program exits.")
        return

    used_clades = args.clades
    if used_clades == ['all']:
        used_clades = clade_names

    if args.force_all_parents and not args.parents.matches(len(used_clades)):
        print("The number of allowed parents, the number of selected clades and the --force-all-parents conflict so that the results must be empty.")
        return

    if args.mutation_threshold < 0.05 or args.mutation_threshold > 1.0 :
        print("mutation-threshold must be between 0.05 and 1.0")
        return

    global reference
    vprint("Reading reference genome, lineage definitions...")
    reference = read_fasta('reference.fasta', None)['MN908947 (Wuhan-Hu-1/2019)']
    all_examples = read_examples('virus_properties.json')

    vprint("Done.\nReading actual input.")
    all_samples = dict()
    for path in args.input:
        read_samples = read_subs_from_fasta(path)
        for key, val in read_samples.items():
            all_samples[key] = val
    vprint("Done.")

    global primer_sets
    primer_sets = dict()
    if args.primers:
        vprint("Reading primers.")
        for path in args.primers:
            pools = read_bed(path)
            primer_sets[path] = pools
        vprint("Done.")

    used_examples = dict()
    for ex_name, ex in all_examples.items():
        if ex_name in used_clades:
            used_examples[ex_name] = ex

    calculate_relations(used_examples)

    match_sets = dict()

    vprint("Scanning input for matches against linege definitons...")
    for sa_name, sa in all_samples.items():
        matching_example_names = []
        if args.force_all_parents:
            matching_example_names = used_examples.keys()
        else:    
            for ex_name, ex in used_examples.items():
                matches_count = len(sa['subs_set'] & ex['unique_subs_set'])
                if matches_count >= args.unique: # theoretically > 0 already gives us recombinants, but they are much more likely to be errors or coincidences
                    matching_example_names.append(ex_name)# 

        matching_examples_tup = tuple(matching_example_names)

        if args.parents.matches(len(matching_example_names)):
            #print(f"{sa_name} is a possible recombinant of {len(matching_example_names)} lineages: {matching_example_names}")
            if match_sets.get(matching_examples_tup):
                match_sets[matching_examples_tup].append(sa)
            else:
                match_sets[matching_examples_tup] = [sa]

    vprint("Done.\nPriniting detailed analysis:\n\n")

    for example_names, samples in match_sets.items():
        show_matches(used_examples, example_names, samples)

def update_readme(parser: argparse.ArgumentParser):
    # on wide monitors, github displays up to 90 columns of preformatted text
    # but 10% of web users have screens which can only fit 65 characters
    global width_override
    width_override = 65
   
    help = parser.format_help()

    new_lines = []

    between_markers = False

    with open("README.md", "rt") as old_readme:
        for line in old_readme:
            if line.strip() == "<!-- BEGIN_MARKER -->":
                between_markers = True
                new_lines.append('<!-- BEGIN_MARKER -->\n')

            if line.strip() == "<!-- END_MARKER -->":
                between_markers = False
                new_lines.append('```\n')
                new_lines.append(help + '\n')
                new_lines.append('```\n')

            if not between_markers:
                new_lines.append(line)

    with open("README.md", "wt") as new_readme:
        new_readme.writelines(new_lines)
            


def vprint(text: str):
    if args.verbose:
        print(text)

def rebuild_examples():
    print("Rebuilding examples from cov-spectrum.org...")
    with open('virus_properties.json', newline='', mode='w') as jsonfile:

        the_map = {}

        for clade, clade_props in mappings['by_clade'].items():
            pango = clade_props['PangoLineage']
            print(f"Fetching data for {clade} / {pango}")
            r = requests.get(f'https://lapis.cov-spectrum.org/open/v1/sample/nuc-mutations?pangoLineage={pango}*&minProportion=0.05')
            result = r.json()
            if len(result['errors']):
                print("Errors occured while querying cov-spectrum.org:")
                for e in result['errors']:
                    print("    " + e)

            the_map[clade] = result['data']

        props = {
            "schemaVersion": "s2r 0.0.1",
            "comment": "File format is inspired by Nextstrains virus_properties, but evolved to be more like cov-spectrum API.",
            "nucMutLabelMap": {},
            'nucMutLabelMapReverse': the_map
        }

        json.dump(props, jsonfile, indent=4)
        print("Examples written to disk.")

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
        assert props['schemaVersion'] == 's2r 0.0.1'
        sequences = {}
        for clade, subs in props['nucMutLabelMapReverse'].items():
            name = clade # TODO get pango-names
            subs_dict = dict()
            for s in subs:
                # I experimented with different file formats, this code can read different ones
                if isinstance(s, str):
                    s = s.strip()
                else:
                    if s['proportion'] < args.mutation_threshold:
                        continue
                    s = s['mutation'].strip()
                
                if len(s) > 0:
                    sub = parse_sub(s)
                    if (sub.mut != '-' or args.enable_deletions) and sub.mut != '.':
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

class Primer(NamedTuple):
    start: int
    end: int
    direction: str
    alt: bool
    name: str
    sequence: str

class Amplicon:
    left_primers: list
    right_primer: list
    number: int
    color: str
    start: int
    end: int

    def __init__(self, number: int):
        self.number = number
        self.left_primers = list()
        self.right_primers = list()
        self.color = get_color(number)
        self.start = None
        self.end = None
        self.amp_start = None
        self.amp_end = None

    def __str__(self):
        return f"Amplicon {self.number} ({self.start} to {self.end})"

    def add_primer(self, primer):
        if primer.direction == '+':
            self.left_primers.append(primer)
            if self.amp_start:
                self.amp_start = max(self.amp_start, primer.end + 1)
            else:
                self.amp_start = primer.end + 1
        else:
            self.right_primers.append(primer)
            if self.amp_end:
                self.amp_end = min(self.amp_end, primer.start - 1)
            else:
                self.amp_end = primer.start - 1
        
        if self.start:
            self.start = min(self.start, primer.start)
        else:
            self.start = primer.start  

        if self.end:
            self.end = max(self.end, primer.end)
        else:
            self.end = primer.end

    def get_char(self, coord: int):
        if coord <= self.start or coord >= self.end:
            return ' '

        for primer in self.left_primers:
            if primer.start <= coord and primer.end >= coord:
                if primer.alt:
                    return '‹'
                else:
                    return '«'

        for primer in self.right_primers:
           if primer.start <= coord and primer.end >= coord:
                if primer.alt:
                    return '›'
                else:
                    return '»'

        return '-'
    
    def overlaps_coord(self, coord: int, actual_amplicon: bool):
        if actual_amplicon:
            return coord >= self.amp_start and coord <= self.amp_end
        else:
            return coord >= self.start and coord <= self.end

    def overlaps_interval(self, interval: Interval):
        if interval.max and self.start > interval.max:
            return False
        if interval.min and self.end < interval.min:
            return False
        return True
        

def read_bed(path):
    pools = dict()
    index = 0
    current_name = None
    with open(path, newline='') as bed:
        for line in bed:
            parts = line.strip().split("\t")

            if(parts[0] == parts[3]): # EasySeq format
                name = parts[6]
                name_parts = name.split("_")
                amplicon_index = int(name_parts[1])
                amplicon = Amplicon(amplicon_index)
                left_primer = Primer(
                    start = int(parts[1]),
                    end = int(parts[2]),
                    name = "left_" + str(amplicon_index),
                    alt = False,
                    direction = "+",
                    sequence = None
                )
                right_primer = Primer(
                    start = int(parts[4]),
                    end = int(parts[5]),
                    name = "right_" + str(amplicon_index),
                    alt = False,
                    direction = "-",
                    sequence = None
                )
                amplicon.add_primer(left_primer)
                amplicon.add_primer(right_primer)

                pool_index = (amplicon_index + 1) % 2 + 1

                if not pools.get(pool_index):
                    pools[pool_index] = dict()

                pools[pool_index][amplicon_index] = amplicon

            else: 
                # ARTIC format
                name = parts[3]
                name_parts = name.split("_")
                amplicon_index = int(name_parts[1])
                pool_index = parts[4]
                direction = parts[5]

                pool = pools.get(pool_index)
                if not pool:
                    pool = dict()
                    pools[pool_index] = pool

                amplicon = pool.get(amplicon_index)
                if not amplicon:
                    amplicon = Amplicon(amplicon_index)
                    pool[amplicon_index] = amplicon

                primer = Primer(
                    start = int(parts[1]),
                    end = int(parts[2]),
                    name = parts[3],
                    alt = len(name_parts) == 4,
                    direction = direction,
                    sequence = parts[6] if 6 < len(parts) else None
                )

                amplicon.add_primer(primer)


    return pools

def read_fasta(path, index_range):
    sequences = dict()
    index = 0
    current_name = None
    with open(path, newline='') as fasta:
        current_sequence = ''
        for line in fasta:
            if(line[0] == '>'):
                if current_name and (not index_range or index_range.matches(index)):
                    sequences[current_name] = current_sequence
                index += 1
                if index_range and index > index_range.max:
                    return sequences
                current_sequence = ''
                current_name = line[1:].strip()
            else:
                current_sequence += line.strip().upper()
        sequences[current_name] = current_sequence

    return sequences

def read_subs_from_fasta(path):
    fastas = read_fasta(path, args.select_sequences)
    sequences = dict()
    start_n = -1
    removed_due_to_ambig = 0
    for name, fasta in fastas.items():
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
 
    if removed_due_to_ambig:
        print(f"Removed {removed_due_to_ambig} of {len(fastas)} sequences with more than { args.max_ambiguous} ambiguous nucs.")

    return sequences

class Sub(NamedTuple):
    ref: str
    coordinate: int
    mut: str

def parse_sub(s):
    if(s[0].isdigit()):
        coordinate = int(s[0:-1])
        return Sub(reference[coordinate-1], coordinate, s[-1])
    else:
        return Sub(s[0], int(s[1:-1]), s[-1])

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

    if args.sort_by_id:
        samples.sort(key = lambda sample: sample['name'][:args.sort_by_id])

    examples = [all_examples[name] for name in example_names]

    coords = set()
    for ex in examples:
        for sub in ex['subs_list']:
            coords.add(sub.coordinate)

    if args.show_private_mutations:
        for sa in samples:
            for sub in sa['subs_list']:
                coords.add(sub.coordinate)

    if args.primers:
        for name, primer_set in primer_sets.items():
            for pool in primer_set.values():
                for amplicon in pool.values():
                    if args.primer_intervals:
                        # check if amplicon should be shown at all, or if it's outside primer_intervals
                        amplicon_matches = False
                        for interval in args.primer_intervals:
                            if amplicon.overlaps_interval(interval):
                                amplicon_matches = True
                                break
                        if not amplicon_matches:
                            continue

                    # check if enough of the actual amplicon range is shown to display its number
                    name_len = len(str(amplicon.number))
                    matched_coords = 0
                    for coord in coords:
                        if amplicon.overlaps_coord(coord, True):
                            matched_coords += 1
                    if matched_coords < name_len:
                        coords.update(range(amplicon.amp_start, amplicon.amp_start + name_len))

                    # make sure that every alt primer is shown for at least one coord
                    # otherwise mismatches in the primary primer may look as if they
                    # would not be compensated by an alt primer
                    for primer in (amplicon.left_primers + amplicon.right_primers):
                        if primer.alt:
                            coords.add(primer.start)

                    # if amplicon.number == 76:
                    #     coords.update(range(amplicon.start, amplicon.end + 1))

    ordered_coords = list(coords)
    ordered_coords.sort()
    
    color_by_name = dict()
    color_index = 0
    for ex in examples:
        color_by_name[ex['name']] = get_color(color_index)
        color_index += 1

    # This method works in a weird way: it pre-constructs the lines for the actual sequences, 
    # and while it constructs the strings, it decides if they are worth showing at the same time.
    # Then, if at least one such string was collected, it prints the header lines for them, and after that the strings.

    ###### SHOW SAMPLES
    current_color = 'grey'
    collected_outputs = []

    last_id = ""

    for sa in samples:
        #current_color = get_color(color_index)
        #color_by_name[sa['name']] = current_color

        prev_definitive_match = None
        breakpoints = 0
        definitives_since_breakpoint = 0
        definitives_count = []

        output = ''

        output += fixed_len(sa['name'], ml) + ' '
        for c, coord in enumerate(ordered_coords):
            if args.add_spaces and c % args.add_spaces == 0:
                output += " "
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
            postfix = '/' + str(num_intermissions)
            num_breakpoints += (num_intermissions - args.max_intermission_count) * 2
            num_intermissions = args.max_intermission_count

        output += f" {num_breakpoints} BP"
        if num_intermissions:
            output += f", {num_intermissions}{postfix} I <= {args.max_intermission_length}"

        if args.breakpoints.matches(num_breakpoints):
            if args.sort_by_id and args.sort_by_id != 999 and last_id != sa['name'][:args.sort_by_id] != last_id[:args.sort_by_id]:
                collected_outputs.append("---")

            last_id = sa['name']
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

            for c, coord in enumerate(ordered_coords):
                if args.add_spaces and c % args.add_spaces == 0:
                    prunt(" ")
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

        for c, coord in enumerate(ordered_coords):
            for name, limits in genes.items():
                if coord >= limits[0] and coord <= limits[1]:
                    if current_name != name:
                        current_name = name
                        color_index += 1
                        current_color = get_color(color_index)
                        text_index = 0
                    
            # Do this once or twice, depending on space insertion
            for i in range(1 + (args.add_spaces and c % args.add_spaces == 0)):
                char = ' '
                if len(current_name) > text_index:
                    char = current_name[text_index]
                cprint(char, 'grey', "on_" + current_color, end='')
                text_index += 1

        print(" ")

        if args.primers:
            ###### SHOW PRIMERS

            prunt('\n')
            for name, primer_set in primer_sets.items():
                for index, pool in primer_set.items():
                    prunt(fixed_len(f"{name}, pool {index}", ml + 1))

                    for c, coord in enumerate(ordered_coords):
                        char = ' '
                        for amplicon in pool.values():

                            if args.primer_intervals:
                                amplicon_matches = False
                                for interval in args.primer_intervals:
                                    if amplicon.overlaps_interval(interval):
                                        amplicon_matches = True
                                        break
                                if not amplicon_matches:
                                    continue

                            if amplicon.overlaps_coord(coord, False):
                                char = amplicon.get_char(coord)
                                if current_name != str(amplicon.number):
                                    current_name = str(amplicon.number)
                                    text_index = 0
                                current_color = amplicon.color

                        if args.add_spaces and c % args.add_spaces == 0:
                            prunt(' ')

                        if char == '-' and len(current_name) > text_index:
                            char = current_name[text_index]
                            text_index += 1
                        cprint(char, current_color, end='')

                    print(' ')
            
                print()

        ###### SHOW REF
        
        prunt(fixed_len("ref", ml + 1))
        for c, coord in enumerate(ordered_coords):
            if args.add_spaces and c % args.add_spaces == 0:
                prunt(" ")
            prunt(reference[coord-1])
        print()
        print()

        ###### SHOW EXAMPLES

        for ex in examples:
            current_color = color_by_name[ex['name']]
            prunt(fixed_len(pretty_name(ex['name']), ml) + ' ', current_color)
            for c, coord in enumerate(ordered_coords):
                if args.add_spaces and c % args.add_spaces == 0:
                    prunt(" ")
                if(ex['subs_dict'].get(coord)):
                    prunt(ex['subs_dict'][coord].mut, current_color)
                else:
                    prunt("•")
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
        vprint(f"Clade  {pretty_name(example['name'])} has {len(example['subs_set'])} mutations, of which {len(example['unique_subs_set'])} are unique.")

class ArgumentAdvancedDefaultsHelpFormatter(argparse.HelpFormatter):
    """In contrast to ArgumentDefaultsHelpFormatter from argparse,
    this formatter also shows 'const' values if they are present, and
    adds blank lines between actions.
    """

    def __init__(self,
                 prog,
                 indent_increment=2,
                 max_help_position=24,
                 width=None):

        global width_override

        if width_override:
            width = width_override

        super().__init__(prog,
                 indent_increment,
                 max_help_position,
                 width)

    def _get_help_string(self, action):
        help = action.help
        if '%(default)' not in action.help and not isinstance(action, argparse._StoreConstAction):
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    if action.const:
                        help += ' (default without flag: %(default)s, default with flag: %(const)s)'
                    else:
                        help += ' (default: %(default)s)'
        return help

    def _format_action(self, action):
        return super()._format_action(action) + '\n'


if __name__ == '__main__':
    main()
