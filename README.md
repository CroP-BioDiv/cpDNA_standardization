# cpDNA standardization

Performs standardization of chloroplast DNA sequence as described in an article
"Towards the Well-Tempered Chloroplast DNA Sequences".

Script is implemented in Python3.

Script converts one or more annotated sequences, in GenBank format.
Converted sequences are in fasta format.

In case of one input sequence, ouput can be specified as specific filename with `-o` argument.
In case of more input sequences, output is specified as directory name, and output files are
named based on input filenames by changing extension into '.fa'

Note: some annotation tools can annotate inverted repeats even if they are not exact repeats,
but support slightly 'error' tolerance. Script's argument `-l` / `--length-ir-difference`
specifies maximal difference in length between IRa and IRb to tolerate. That means if length
argument is set, annotations where length of repeats differ more then specified number will be
discarded, and sequence will be threated as IRs are not annotated at all. Article
"Towards the Well-Tempered Chloroplast DNA Sequences" for this argument uses value 10 (bp).

Usage:
```
# Help
python3 cpdna_standardization.py -h

# Convert one sequence
python3 cpdna_standardization.py <path_to_gen_bank_seq> -o output_sequence.fa

# Convert more sequences
python3 cpdna_standardization.py <dir_to_gen_bank_seqs>/*.gb -d output_sequences
```
