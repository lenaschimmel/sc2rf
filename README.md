# What's this?
This program can search genome sequences of SARS-CoV-2 for potential recombinants - new virus lineages that have (partial) genes from more than one parent lineage.

**The initial version of this program has been written as a quick-and-dirty proof-of-concept and needs a lot of cleanup.**

# Requirements and Installation
You need at least Python 3.6 and you need to install the requirements first. You might use somethink like `pip3 install -r requirements.txt` to do that.

Also, you need a terminal which supports ANSI control sequences. On Linux, MacOS, etc. it should probably work. On Windows, you need a recent version of Windows 10 and you run the script from `cmd.exe`. See [this table](https://pypi.org/project/termcolor/) for details.

# Usage
Start with a `sample.fasta` file with one or more sequences which might contain recombinants. (Of course, future versions will allow you to name your input file however you like, etc.) Your sequences have to be aligned to the `reference.fasta`. If they are not, you will get an error message like:

> Sequence hCoV-19/Phantasialnd/EFWEFWD not properly aligned, length is 29718 instead of 29903.

_(For historical reasons, I always used [Nextclade](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli.html) to get aligned sequences, but you might also use [Nextalign](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextalign-cli.html) or any other tool.)_

Run this script, like:

```bash
python3 search_recombinants.py
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
   * [ ] as input file
   * [ ] as piped stream
 * [ ] If we still accept csv/ssv input, autodetect the delimiter either by file name or by analysing the first line
 * [ ] get rid of `termcolor` which is out only dependency right now
 * [ ] add new pango lineages that do not have their own clade, especially BA.3
 * [ ] find a way to handle already designated recombinant lineages
 * [ ] accept command line arguments to control how it works
   * [ ] min_defining_matches
   * [ ] min_breakpoints
   * [ ] max_breakpoints
   * [ ] min_streak_length
