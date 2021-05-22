import os
from cpdna_standardization import standardize_file, SeqIO


def test_project(project, stop_on_error, length_ir_difference):
    if not os.path.isdir(project):
        print(f"Project directory ({project}) doesn't exist!")
        return
    in_dir = os.path.join(project, '03_GeSeq')
    if not os.path.isdir(in_dir):
        print(f"Input directory ({in_dir}) doesn't exist!")
        return
    out_dir = os.path.join(project, 'Gn_01_seq')
    if not os.path.isdir(out_dir):
        print(f"Output directory ({out_dir}) doesn't exist!")
        return

    errors = []
    for f_gb in os.listdir(in_dir):
        if f_gb.endswith('.gb'):
            output_seq_rec = standardize_file(os.path.join(in_dir, f_gb), length_ir_difference=length_ir_difference)
            #
            p_file = os.path.join(out_dir, f"p_{f_gb[3:-3]}.fa")
            gb_file = os.path.join(out_dir, f"p_{f_gb[3:]}")
            gb_exists = os.path.isfile(gb_file)
            fa_exists = os.path.isfile(p_file)

            if output_seq_rec is None:
                if gb_exists or fa_exists:
                    errors.append(f"Error: no annotation, but file in output directory ({f_gb})!")
                    if stop_on_error:
                        break

            elif output_seq_rec is True:
                if not os.path.join(out_dir, f_gb):
                    errors.append(f"Error: standardized, but gb file not in output directory ({f_gb})!")
                    if stop_on_error:
                        break
                if gb_exists or fa_exists:
                    errors.append(f"Error: standardized, but p_ file in output directory ({f_gb})!")
                    if stop_on_error:
                        break

            else:
                if not gb_exists and not fa_exists:
                    errors.append(f"Error: not standardized, but no file in output directory ({f_gb})!")
                    if stop_on_error:
                        break
                else:
                    st_seq_rec = SeqIO.read(gb_file, 'genbank') if gb_exists else SeqIO.read(p_file, 'fasta')
                    if str(output_seq_rec.seq) != str(st_seq_rec.seq):
                        errors.append(f"Error: not standardized, fa differs ({f_gb})!")
                        if stop_on_error:
                            break
    #
    return errors


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="""
Tests chloroplast project data.
Input directory with annotated sequences is subdirectory 03_GeSeq. Filenames are names NC_<num>.gb.
Standardized sequences are in subdirectory Gn_01_seq. Filenames are names p_<num>.fa, p_<num>.gb, and NC_<num>.gb.
""")
    parser.add_argument('project', help='Project directory.')
    parser.add_argument('-s', '--stop-on-error', action='store_true', help='Stop processing on error')
    parser.add_argument('-l', '--length-ir-difference', default=10, type=int,
                        help='Max difference in length between IRa and IRb')
    params = parser.parse_args()

    if (errors := test_project(params.project, params.stop_on_error, params.length_ir_difference)):
        for e in errors:
            print(e)
