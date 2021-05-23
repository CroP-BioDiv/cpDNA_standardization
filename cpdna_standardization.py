import os
from Bio import SeqIO, Seq, SeqRecord
from Bio.SeqFeature import FeatureLocation, CompoundLocation


def find_chloroplast_irs(seq_rec, length_ir_difference=None):
    # Finds the longest pair of inverted repeats
    # Note: in 
    _ir = ('inverted',)
    rep_regs = [f for f in seq_rec.features
                if f.type == 'repeat_region' and f.qualifiers.get('rpt_type', _ir)[0] == 'inverted']

    if len(rep_regs) < 2:
        return

    # Repair features with location of type <~start..ira.start>
    len_seq = len(seq_rec)
    half_length = len_seq // 2
    for f in rep_regs:
        if len(f) >= half_length:
            assert False, seq_rec.id  # For test
            loc = f.location
            f.location = CompoundLocation([FeatureLocation(loc.end, len_seq, strand=1),
                                           FeatureLocation(0, loc.start + 1, strand=1)])

    if length_ir_difference is not None:
        max_len = max(map(len, rep_regs)) - length_ir_difference
        max_regs = [f for f in rep_regs if len(f) >= max_len]
        if len(max_regs) != 2:
            return None
        ira, irb = max_regs
    else:
        ira, irb = sorted(rep_regs, key=len)[-2:]

    # Take two the longest regions
    ira, irb = sorted(rep_regs, key=len)[-2:]
    diff_1 = (irb.location.parts[0].start - ira.location.parts[-1].end) % len_seq
    diff_2 = (ira.location.parts[0].start - irb.location.parts[-1].end) % len_seq
    return (ira, irb) if diff_1 < diff_2 else (irb, ira)


def create_chloroplast_partition(len_seq, ira_s, ira_e, irb_s, irb_e):
    # Partition is a list of tuples (region_name, from_idx, to_idx).
    # It is possible that last and first region_name are the same
    partition = []

    def _add_ir(name, loc_start, loc_end):
        if loc_end > loc_start:
            partition.append((name, loc_start, loc_end))
        else:
            partition.append((name, loc_start, len_seq))
            partition.append((name, 0, loc_end))

    _add_ir('ira', ira_s, ira_e)
    _add_ir('ssc', ira_e, irb_s)
    _add_ir('irb', irb_s, irb_e)
    _add_ir('lsc', irb_e, ira_s)
    return partition


def _intersect_intervals(s1, e1, s2, e2):
    return (s1 <= s2 < e1) or \
        (s1 < e2 <= e1) or \
        (s2 <= s1 < e2)


def _intersect(feature, from_idx, to_idx):
    return any(_intersect_intervals(from_idx, to_idx, p.start, p.end) for p in feature.location.parts)


def _copy_sequence_annotations(from_seq_rec, to_seq_rec):
    # Fixes: "ValueError: missing molecule_type in annotations"
    f_ann = from_seq_rec.annotations
    t_ann = to_seq_rec.annotations
    for k in ('molecule_type', 'topology', 'accessions', 'organism', 'taxonomy'):
        if v := f_ann.get(k):
            t_ann[k] = v


def chloroplast_parts_orientation(seq_rec, partition):
    # Check chloroplast sequence part orientation.
    # Default orientation is same as one uses in Fast-Plast. Check:
    #  - source file orientate_plastome_v.2.0.pl
    #    (https://github.com/mrmckain/Fast-Plast/blob/master/bin/orientate_plastome_v.2.0.pl)
    #  - explanation https://github.com/mrmckain/Fast-Plast/issues/22#issuecomment-389302606
    # Consitent with Wikipedia image:
    #  - https://en.wikipedia.org/wiki/Chloroplast_DNA#/media/File:Plastomap_of_Arabidopsis_thaliana.svg

    l_seq = len(seq_rec)
    in_parts = dict((r, []) for r in ('lsc', 'ira', 'ssc', 'irb'))  # lists of tuples (feature name, feature strand)
    for f in seq_rec.features:
        if f.type == 'gene' and f.location:
            in_regs = set(reg_name for reg_name, from_idx, to_idx in partition if _intersect(f, from_idx, to_idx))
            f_name = f.qualifiers['gene'][0]
            if not in_regs:
                print(f'Warning: gene {f_name} is not in any region!')
            elif len(in_regs) > 2:
                print(f'Warning: gene {f_name} in more than two regions!', in_regs)
            else:
                for r in in_regs:
                    in_parts[r].append((f_name, f.strand))

    lsc_count = sum(s for n, s in in_parts['lsc'] if any(x in n for x in ('rpl', 'rps')) and s)
    ssc_count = sum(s for n, s in in_parts['ssc'] if s)
    ira_count = sum(s for n, s in in_parts['ira'] if 'rrn' in n and s)

    return dict(lsc=(lsc_count <= 0),
                ssc=(ssc_count <= 0),
                ira=(ira_count >= 0))


def standardize(seq_rec, length_ir_difference=None, info=False):
    # Returns:
    #  - None if it IR annotation is missing,
    #  - True if sequence is in standard form
    #  - SeqRecord object of seqeunce in standard form
    if not (irs := find_chloroplast_irs(seq_rec, length_ir_difference=length_ir_difference)):
        if info:
            print(f"info: sequence {seq_rec.id} doesn't have annotated IRs!")
        return

    len_seq = len(seq_rec)
    ira, irb = irs
    ira_s, ira_e = ira.location.parts[0].start, ira.location.parts[-1].end
    irb_s, irb_e = irb.location.parts[0].start, irb.location.parts[-1].end
    partition = create_chloroplast_partition(len_seq, ira_s, ira_e, irb_s, irb_e)
    region_orientations = chloroplast_parts_orientation(seq_rec, partition)

    seq_rec.features = []  # Remove features, since we are storing fasta file

    if all(region_orientations.values()):
        # Regions are good oriented, check cyclic shift
        if irb_e <= 10 or (len_seq - irb_e) <= 10:
            if info:
                print(f"info: sequence {seq_rec.id} is in standard form.")
            return True
        #
        if info:
            print(f'info: sequence {seq_rec.id}, cycle shifted for {irb_e}. Sequence length {len_seq}.')
        new_seq_rec = seq_rec[irb_e:] + seq_rec[:irb_e]
    else:
        seq = str(seq_rec.seq)
        regions = dict((n, '') for n in ('lsc', 'ira', 'ssc', 'irb'))
        for n, f, t in partition:
            regions[n] += seq[f:t]
        regions = dict((n, SeqRecord.SeqRecord(Seq.Seq(dna))) for n, dna in regions.items())

        if not region_orientations['lsc']:  # LSC
            regions['lsc'] = regions['lsc'].reverse_complement()
        if not region_orientations['ssc']:  # SSC
            regions['ssc'] = regions['ssc'].reverse_complement()
        if not region_orientations['ira']:  # IRs
            regions['ira'], regions['irb'] = regions['irb'].reverse_complement(), regions['ira'].reverse_complement()

        if info:
            o_regions = ", ".join(n for n, o in region_orientations.items() if not o)
            print(f'info: sequence {seq_rec.id}, reoriented regions: {o_regions}')
        new_seq_rec = regions['lsc'] + regions['ira'] + regions['ssc'] + regions['irb']

    #
    _copy_sequence_annotations(seq_rec, new_seq_rec)
    return new_seq_rec


def standardize_file(input_filename, length_ir_difference=None, info=False):
    input_seq_rec = SeqIO.read(input_filename, 'genbank')
    return standardize(input_seq_rec, length_ir_difference=length_ir_difference, info=info)


def _convert(input_filename, output_filename, length_ir_difference):
    output_seq_rec = standardize_file(input_filename, length_ir_difference=length_ir_difference, info=True)
    if output_seq_rec and output_seq_rec is not True:
        SeqIO.write([output_seq_rec], open(output_filename, 'w'), 'fasta')


def run(input_filename, output_filename, output_dirname, length_ir_difference):
    if len(input_filename) == 1 and output_filename:
        _convert(input_filename[0], output_filename, length_ir_difference)
    else:
        if not os.path.isdir(output_dirname):
            os.makedirs(output_dirname)

        for f in input_filename:
            base_name = os.path.splitext(os.path.basename(f))[0]
            _convert(f, os.path.join(output_dirname, f'{base_name}.fa'), length_ir_difference)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="""
Script converts one or more annotated sequences, in GenBank format.
Converted sequences are in fasta format.

In case of one input sequence, ouput can be specified as specific filename with `-o` argument.

In case of more input sequences, output is specified as directory name, and output files are
named based on input filenames by changing extension into '.fa'
""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('input_filename', nargs='+', help='Input filename(s).')
    parser.add_argument('-o', '--output-filename', help='Output filename.')
    parser.add_argument('-d', '--output-dirname', default='.',
                        help='Output directory. Default current working directory.')
    parser.add_argument('-l', '--length-ir-difference', default=None, type=int,
                        help='Max difference in length between IRa and IRb.')
    # ToDo:
    # parser.add_argument('-s', '--save-standardized', help='Save standardized also')
    # parser.add_argument('-n', '--save-not-annotated', help='Save not annotated')
    params = parser.parse_args()

    run(params.input_filename, params.output_filename, params.output_dirname, params.length_ir_difference)
