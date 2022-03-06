# What's this?
This program can search genome sequences of SARS-CoV-2 for potential recombinants - new virus lineages that have (partial) genes from more than one parent lineage.

**The initial version of this program has been written as a quick-and-dirty proof-of-concept and needs a lot of cleanup.**

# Requirements and Installation
You need at least Python 3.6 and you need to install the requirements first. You might use somethink like `pip3 install -r requirements.txt` to do that.

# Usage
Start with a `.fasta` file with one or more sequences which might contain recombinants.

Currently, you need to run your sequences through Nextclade before sarscov2recombinants can use them. This preparation step will probably become unnecessary very soon, so that you will not need the Nextclade software at all.

```bash
nextclade --in-order --input-fasta sample.fasta --input-dataset path/to/your/datasets/sars-cov-2  --output-dir nextclade-output/  --output-csv sample.ssv --verbose
```

The only file we need is `sample.ssv`. Note how the file is name `.ssv` because it is separated by semicolons, although the option is called `--output-csv`. You can delete everything inside `nextclade-output/` if you want.

Then run this script, like this:

```bash
python3 search_recombinants.py
```

# Source material
I used the following files from Nextclade / Nextstrain:
 * `virus_properties.json`
 * `reference.fasta`

The file `mapping.csv` is a modified version of the table on the [covariants homepage](https://covariants.org/) by Nextstrain.

The initial version of this program was written in cooperation with [@flauschzelle](https://github.com/flauschzelle).

# TODO
 * [ ] provide a sample file (maybe both `.fasta` and `.csv`, as long as the csv step is still needed)
 * [X] don't use sequences as examples - use the lineage defintions from a json file instead
 * [ ] accept aligned fasta 
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