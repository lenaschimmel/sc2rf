# What's this?
This program can search genome sequences of SARS-CoV-2 for potential recombinants - new virus lineages that have (partial) genes from more than one parent lineage.

**The initial version of this program has been written as a quick-and-dirty proof-of-concept and needs a lot of cleanup.**

# Requirements and Installation
You need at least Python 3.6 and you need to install the requirements first. You might use somethink like `pip3 install -r requirements.txt` to do that.

Also, you need a terminal which supports ANSI control sequences. On Linux, MacOS, etc. it should probably work. On Windows, you need a recent version of Windows 10 and you run the script from `cmd.exe`. See [this table](https://pypi.org/project/termcolor/) for details.

# Basic Usage
Start with a `.fasta` file with one or more sequences which might contain recombinants. Your sequences have to be aligned to the `reference.fasta`. If they are not, you will get an error message like:

> Sequence hCoV-19/Phantasialnd/EFWEFWD not properly aligned, length is 29718 instead of 29903.

_(For historical reasons, I always used [Nextclade](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli.html) to get aligned sequences, but you might also use [Nextalign](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextalign-cli.html) or any other tool.)_

Then call:

```
search_recombinants.py <your_filename.fasta>
```

# Advanced Usage
You can execute `search_recombinants.py -h` to get a help message like this one, probably even more up-to-date:

```
usage: search_recombinants.py [-h] [--parents INTERVAL] [--breakpoints INTERVAL]
                              [--clades [{all,20I,20H,20J,21A,21I,21J,21B,21C,21D,...} ...]]
                              [--unique NUM] [--max-intermission-length NUM] [--max-intermission-count NUM]
                              [--max_name_length NUM]
                              input [input ...]

Analyse SARS-CoV-2 sequences for potential, unknown recombinant variants.

positional arguments:
  input                 input sequences to test, as aligned .fasta file(s)

optional arguments:
  -h, --help            show this help message and exit
  --parents INTERVAL    Allowed umber of potential parents of a recombinant. Interval (see below).
  --breakpoints INTERVAL
                        Allowed number of breakpoints in a recombinant. Interval (see below).
  --clades [{all,20I,20H,20J,21A,21I,21J,21B,21C,21D,...} ...]
                        List of clades which are considered as potential parents. Use Nextclade names, i.e. "21A".
                        Also accepts "all".
  --unique NUM          Minimum of substitutions in a sample which are unique to a potential parent clade, so that
                        the clade will be considered.
  --max-intermission-length NUM, -l NUM
                        The maximum length of an intermission in consecutive substitutions. Intermissions are
                        stretches to be ignored when counting breakpoints.
  --max-intermission-count NUM, -m NUM
                        The maximum number of intermissions which will be ignored. Surplus intermissions count
                        towards the number of breakpoints.
  --max_name_length NUM, -n NUM
                        Only show up to NUM characters of sample names.

An Interval can be a single number ("3"), a closed interval ("2-5" ) or an open one ("4-" or "-7"). The limts are
inclusive. Only positive numbers are supported.
```

# Interpreting the output
_To be written..._

Some example output (based on Sequences published by the [German Robert-Koch-Institut](https://github.com/robert-koch-institut/SARS-CoV-2-Sequenzdaten_aus_Deutschland)):

<img width="1110" alt="Screenshot of the terminal output of this program" src="https://user-images.githubusercontent.com/1325019/156946733-cdc025d7-869a-4ce6-b1b7-62b0d1a30bac.png">


# Source material
I used the following files from Nextstrain's [nextclade_data](https://github.com/nextstrain/nextclade_data/tree/master/data/datasets/sars-cov-2/references/MN908947/versions/2022-03-04T12:00:00Z/files):
 * `virus_properties.json`
 * `reference.fasta`

The file `mapping.csv` is a modified version of the table on the [covariants homepage](https://covariants.org/) by Nextstrain.

The initial version of this program was written in cooperation with [@flauschzelle](https://github.com/flauschzelle).

# TODO
 * [ ] provide a sample file (maybe both `.fasta` and `.csv`, as long as the csv step is still needed)
 * [X] don't use sequences as examples - use the lineage defintions from a json file instead
 * [X] accept aligned fasta 
   * [x] as input file
   * [ ] as piped stream
 * [ ] If we still accept csv/ssv input, autodetect the delimiter either by file name or by analysing the first line
 * [ ] get rid of `termcolor` which is out only dependency right now
 * [ ] add new pango lineages that do not have their own clade, especially BA.3
 * [ ] find a way to handle already designated recombinant lineages
 * [x] accept command line arguments to control how it works
   * [x] min_defining_matches
   * [x] min_breakpoints
   * [x] max_breakpoints
   * [x] min_streak_length
