"""Compares VCF files obtained from a parent-offspring trio.
Produces statistics and variant comparison files, with analysis
of inheritance."""

import os
import sys
import argparse
import matplotlib.pyplot as plt
import vcf


class VcfFile(object):
    """Contains info about VCF file, and current position in file."""
    def __init__(
        self, file_name=None, file_=None, counts=None,
        reader=None, curr_rec=None, next_rec=None, curr_call=None,
        eof=False, has_curr_variant='absent', family_rel=None,
        curr_qc_values=None
    ):
        self.file_name = file_name
        self.file_ = file_
        self.counts = counts    # Category:count map
        self.reader = reader    # PyVCF reader
        self.curr_rec = curr_rec    # Current record
        self.next_rec = next_rec    # Next record
        self.curr_call = curr_call   # Current call
        self.eof = eof    # End of file has been reached
        self.has_curr_variant = has_curr_variant   # Absent/low_qual/high_qual
        self.family_rel = family_rel  # Parent 1/Parent 2/Child
        self.curr_qc_values = curr_qc_values  # DP, QUAL, etc. for current var

    def __str__(self):
        return (
            'Filename: %s\nCounts: %s\nCurr rec: %s\nNext rec: %s\nEOF: %s'
            '\nHas curr variant: %s\nCurr QC values: %s\n' %
            (
                self.file_name, self.counts, self.curr_rec, self.next_rec,
                self.eof, self.has_curr_variant, self.curr_qc_values
            )
        )


def main():
    parser = argparse.ArgumentParser()
    args = parse_arguments(parser)
    validate(args, parser)
    if (
        args.parent_1_file_names is None or
        args.parent_2_file_names is None or
        args.child_file_names is None
    ):
        raise ValueError(
            'Parent 1, Parent 2, and Child file names must be specified.'
        )

    # Per-file statistics, in the order in which they will be displayed
    # in HTML table. (Categories starting with a capital letter are
    # displayed.) Stored in counts dict in VcfFile object.
    counts_per_file_categs = [
        'Total variants', 'Average DP', 'Average QUAL',
        'High-quality variants', 'Homozygous reference', 'SNP',
        'Heterozygous', 'Homozygous alternate', 'Indel', 'Missing',
        'A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T', 'G>A', 'G>C', 'G>T', 'T>A',
        'T>C', 'T>G', 'Predicted sex', 'total_x_dp', 'total_y_dp',
        'total_x_gt', 'total_y_gt', 'count_x', 'count_y',
    ]

    # Categories displayed in flowchart in HTML file
    overall_counts_categs = [
        'Pass QC', 'Indel', 'No indel', 'Internal conflict', 'Invalid',
        'De novo', 'Inherited',
        'Parents HOM-REF, HET; child HOM-REF',
        'Parents HOM-REF, HET; child HET',
        'Parents HOM-REF, HOM-ALT; child HET',
        'Parents HET, HET; child HOM-REF',
        'Parents HET, HET; child HET',
        'Parents HET, HET; child HOM-ALT',
        'Parents HET, HOM-ALT; child HET',
        'Parents HET, HOM-ALT; child HOM-ALT',
        'Parents HOM-ALT, HOM-ALT; child HOM-ALT'
    ]
    overall_counts = {}

    # Categories displayed in QC plots
    qc_counts_categs = [
        'All', 'De novo', 'Inherited', 'Invalid', 'Internal conflict'
    ]
    # QC plot x-axis intervals
    qc_intervals = {
        'dp': xrange(0, 55, 5), 'qd': xrange(0, 20, 2),
        'qual': xrange(0, 260, 26)
    }
    qc_counts = {}

    # Initialize readers, counts
    multi_sample_files = []  # Names of multi-sample VCF files
    vcf_file_objs = []
    for f_name in (
        args.parent_1_file_names + args.parent_2_file_names +
        args.child_file_names
    ):
        in_f = open(f_name)
        reader = vcf.Reader(in_f)
        f_base_name = os.path.basename(f_name)
        if len(reader.samples) > 1:
            multi_sample_files.append(f_base_name)
            in_f.close()
        else:
            v = VcfFile(
                file_=in_f, counts={}, reader=reader, next_rec=reader.next(),
                curr_qc_values={}
            )
            if f_name in args.parent_1_file_names:
                v.family_rel = 'Parent 1'
            elif f_name in args.parent_2_file_names:
                v.family_rel = 'Parent 2'
            else:
                v.family_rel = 'Child'
            for cat in counts_per_file_categs:
                v.counts[cat] = 0
            v.file_name = '%s (%s)' % (f_base_name, v.family_rel)
            vcf_file_objs.append(v)

    if args.create_qc_plots:
        for i in qc_intervals:
            qc_counts[i] = {}
            for j in qc_counts_categs:
                qc_counts[i][j] = {}
                for k in qc_intervals[i]:
                    qc_counts[i][j][k] = 0
    for x in overall_counts_categs:
        overall_counts[x] = 0

    # Iterate over VCF files, write variant comparison output files,
    # fill in counts
    header = (
        'CHROM\tPOS\tREF\tALT\tGT\tInheritance_category\t'
        'Inheritance_subcategory\t%s\t%%Present\tAvg_QD'
        '\tAvg_DP\tAvg_QUAL\tAvg_AD0\tAvg_AD1\n' %
        ('\t'.join([v.file_name for v in vcf_file_objs]))
    )
    with open(args.out_all_variants, 'w') as out_all:
        out_all.write(header)
        compare_variants(
            vcf_file_objs, qc_counts, overall_counts, args, out_all
        )

    for v in vcf_file_objs:
        v.file_.close()
        for category in v.counts:
            if 'average' in category.lower():
                v.counts[category] /= float(v.counts['Total variants'])
        v.counts['Predicted sex'] = predict_sex(
            v.counts['total_x_dp'], v.counts['total_y_dp'],
            v.counts['total_x_gt'], v.counts['total_y_gt'],
            v.counts['count_x'], v.counts['count_y']
        )

    # Output QC plots
    if args.create_qc_plots:
        output_qc_plot(
            qc_counts_categs, 'Read depth', qc_counts['dp'], args.depth_plot,
            qc_intervals['dp']
        )
        output_qc_plot(
            qc_counts_categs, 'Quality by depth', qc_counts['qd'],
            args.qd_plot, qc_intervals['qd']
        )
        output_qc_plot(
            qc_counts_categs, 'Quality', qc_counts['qual'], args.qual_plot,
            qc_intervals['qual']
        )

    # Output HTML stats file
    counts_norm = normalize_counts(
        counts_per_file_categs, [v.counts for v in vcf_file_objs]
    )
    with open(args.out_stats, 'w') as out_f:
        page_title = 'VCF Comparison'
        out_f.write(
            '<!DOCTYPE html><html><head><title>%s</title></head><body>'
            '<h1>%s</h1>' % (page_title, page_title)
        )
        table_message = (
            '<p>Counts per unique combination of CHROM and POS. </p>'
            '<p>Normalized with respect to the number of '
            'high-quality variants, and to the value in the leftmost '
            'numeric column. Normalized values less than 0.95 or '
            'greater than 1.05 are shown in blue. Counts shown '
            '(excluding total variants) are after quality filtering. '
            '</p>'
        )

        print_dict_list(
            out_f, 'Counts Per File', counts_norm,
            [z for z in counts_per_file_categs if z[0].isupper()],
            [v.file_name for v in vcf_file_objs], table_message
        )
        if len(multi_sample_files) > 0:
            out_f.write(
                '<p><b>Multi-sample files excluded from '
                'analysis: </b>%s</p>' % (', '.join(multi_sample_files))
            )
        out_f.write(
            '<h2>Overall Counts Per Unique Combination of CHROM, POS, '
            'REF, ALT, and GT</h2>'
        )
        out_f.write(
            '<p>Each category (except Pass QC) is a subcategory of the '
            'rightmost category directly above it.</p>'
        )
        out_f.write(
            '<table style="table-layout:fixed; width:'
            '%dpx; word-wrap:break-word"><tr>' % 2000
        )
        for x in overall_counts_categs:
            if x in [
                'Indel', 'Internal conflict',
                'Parents HOM-REF, HET; child HOM-REF'
            ]:
                out_f.write('</tr><tr>')
            out_f.write(
                '<td style=padding-right:2em><b>%s</b>: %d</td>' %
                (x, overall_counts[x])
            )
        out_f.write('</tr></table>')
        out_f.write('</body></html>')


def parse_arguments(parser):
    """Adds arguments to parser. Returns parsed arguments."""
    parser.add_argument(
        '-p1', '--parent_1_file_names', nargs='+',
        help='VCF filenames corresponding to Parent 1 (required)'
    )
    parser.add_argument(
        '-p2', '--parent_2_file_names', nargs='+',
        help='VCF filenames corresponding to Parent 2 (required)'
    )
    parser.add_argument(
        '-c', '--child_file_names', nargs='+',
        help='VCF filenames corresponding to child (required)'
    )
    parser.add_argument(
        '--create_qc_plots', action='store_true',
        help='Whether to output QC plots. Takes no arguments. '
             '(Default: Does not output plots.)'
    )
    parser.add_argument(
        '--out_stats', default='vcf_stats.html',
        help='Name of HTML table output file containing overall '
             'statistics. (Default: vcf_stats.html)'
    )
    parser.add_argument(
        '--depth_plot', default='dp_plot.png',
        help='Name of output file containing read depth QC plot. '
             '(--create_qc_plots must be specified.) (Default: '
             'dp_plot.png)'
    )
    parser.add_argument(
        '--qd_plot', default='qd_plot.png',
        help='Name of output file containing quality by depth QC '
             'plot. (--create_qc_plots must be specified.) (Default: '
             'qd_plot.png)'
    )
    parser.add_argument(
        '--qual_plot', default='qual_plot.png',
        help='Name of output file containing quality QC plot. '
             '(--create_qc_plots must be specified.) (Default: '
             'qual_plot.png)'
    )
    parser.add_argument(
        '--out_all_variants', default='all_variants.txt',
        help='Name of output file containing all high-quality '
             'variants. (Default: all_variants.txt)'
    )
    parser.add_argument(
        '--min_qd', type=int, default=3,
        help='Minimum quality by depth for inclusion in variant '
        'output files. (Default: 3)'
    )
    parser.add_argument(
        '--min_qual', type=int, default=80,
        help='Minimum quality score for inclusion in variant output '
             'files. (Default: 80)'
    )
    parser.add_argument(
        '--min_dp', type=int, default=15,
        help='Minimum read depth for inclusion in variant output '
             'files. (Default: 15)'
    )
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()


def validate(args, parser):
    """Validate command line arguments."""
    message = ''
    for x in (
        args.parent_1_file_names + args.parent_2_file_names +
        args.child_file_names
    ):
        message += check_extension(x, '.vcf')
    message += check_extension(args.out_stats, '.html')
    message += check_extension(args.out_all_variants, '.txt')

    if message != '':  # Error has occurred
        parser.print_help()
        raise ValueError(message[:-1])


def check_extension(filename, extension):
    """Return error string if wrong extension."""
    if not filename.endswith(extension):
        return 'File "%s" does not have %s extension. ' % (
            filename, extension
        )
    else:
        return ''


def compare_variants(
    vcf_file_objs, qc_counts, overall_counts, args, out_all
):
    """Iterate through all VCF files simultaneously. Write per-file
    variant info to output files and update overall statistics
    (counts).
    """
    for v in vcf_file_objs:
        if (v.reader is None) or (v.next_rec is None) or (v.counts == {}):
            return

    min_chrom = min_chromosome(
        [v.next_rec.CHROM for v in vcf_file_objs if not v.eof]
    )
    print 'Processing %s...' % min_chrom
    while not (all([v.eof for v in vcf_file_objs])):
        new_min_chrom = min_chromosome(
            [v.next_rec.CHROM for v in vcf_file_objs if not v.eof]
        )
        if new_min_chrom != min_chrom:
            print 'Processing %s...' % new_min_chrom
            min_chrom = new_min_chrom
        min_pos = min(
            [
                v.next_rec.POS for v in vcf_file_objs if
                v.next_rec.CHROM == min_chrom and not v.eof
            ]
        )
        for v in vcf_file_objs:
            v.curr_rec = v.next_rec
            v.curr_call = None
            record = v.next_rec
            call = record.genotype(v.reader.samples[0])
            if (
                record.CHROM == min_chrom and record.POS == min_pos
                and not v.eof
            ):
                assign_qc_values(v.curr_qc_values, record, call)
                update_counts_all(v.curr_qc_values, v.counts, record, call)
                if (
                    v.curr_qc_values['qd'] >= args.min_qd and
                    v.curr_qc_values['dp'] >= args.min_dp and
                    v.curr_qc_values['qual'] >= args.min_qual
                ):
                    v.has_curr_variant = 'high_qual'
                    update_counts_high_qual(v.counts, record, call)
                else:
                    v.has_curr_variant = 'low_qual'
                v.curr_call = call
                try:
                    v.next_rec = v.reader.next()
                except StopIteration:
                    v.eof = True
            else:
                v.has_curr_variant = 'absent'

        # Write single variant (at min_chrom and min_pos) to output files,
        # update counts
        process_variant(
            vcf_file_objs, qc_counts, overall_counts, args, out_all
        )


def assign_qc_values(qc_values_dict, record, call):
    """Assign QC values to dict, given current record and call."""
    qc_values_dict['qd'] = record.INFO['QD']
    qc_values_dict['dp'] = call['DP']
    qc_values_dict['qual'] = record.QUAL
    qc_values_dict['ad0'] = call['AD'][0]
    qc_values_dict['ad1'] = call['AD'][1]


def update_counts_all(qc_values, counts, record, call):
    """Update counts dict for a variant."""
    counts['Total variants'] += 1
    counts['Average DP'] += qc_values['dp']
    counts['Average QUAL'] += qc_values['qual']
    if 'X' in record.CHROM.upper():
        counts['total_x_dp'] += qc_values['dp']
        counts['total_x_gt'] += (0 if (call.gt_type == 1) else 1)
        counts['count_x'] += 1
    elif 'Y' in record.CHROM.upper():
        counts['total_y_dp'] += qc_values['dp']
        counts['total_y_gt'] += (0 if (call.gt_type == 1) else 1)
        counts['count_y'] += 1


def update_counts_high_qual(counts, record, call):
    """Update counts dict for a high-quality variant."""
    counts['High-quality variants'] += 1
    if call.gt_type == 1:
        counts['Heterozygous'] += 1
    elif call.gt_type == 2:
        counts['Homozygous alternate'] += 1
    elif call.gt_type is None:
        counts['Missing'] += 1
    else:
        counts['Homozygous reference'] += 1

    counts['Indel'] += record.is_indel
    is_snp = (1 if (record.is_snp and call.is_variant) else 0)
    counts['SNP'] += is_snp
    if is_snp:
        counts[
            '%s>%s' % (
                record.alleles[0],
                record.alleles[int(max(call.gt_alleles))]
            )
        ] += 1


def process_variant(
    vcf_file_objs, qc_counts, overall_counts, args, out_all
):
    """Process a variant uniquely identified by CHROM and POS.
    Find inheritance. Identify and process sub-variants."""
    high_qual_present = any(
        v.has_curr_variant == 'high_qual' for v in vcf_file_objs
    )
    files_per_variant = {}     # (REF, ALT, GT): list of VCF file objs
    if high_qual_present:
        for v in [z for z in vcf_file_objs if z.has_curr_variant != 'absent']:
            call = v.curr_call
            record = v.curr_rec
            ref = str(record.REF)
            alt = str(record.ALT[0])
            gt = call['GT']
            if (ref, alt, gt) in files_per_variant:
                files_per_variant[(ref, alt, gt)].append(v)
            else:
                files_per_variant[(ref, alt, gt)] = [v]

    # Find internal disagreements within Parent 1 files/
    # Parent 2 files/ Child files
    family_rels = ['Parent 1', 'Parent 2', 'Child']
    example_files = {}
    conflict_exists = False
    for family_rel in family_rels:
        example_files[family_rel] = None
        for x in files_per_variant:
            for v in files_per_variant[x]:
                if (
                    v.family_rel == family_rel and
                    v.has_curr_variant != 'absent'
                ):
                    example_files[family_rel] = v
        if internal_conflict_exists(
            vcf_file_objs, family_rel, files_per_variant
        ):
            conflict_exists = True

    # Find inheritance category and subcategory
    subcat = 'NA'
    if True in [
        v.curr_rec.is_indel for v in vcf_file_objs if
        v.has_curr_variant != 'absent'
    ]:
        cat = '[Indel]'
    elif conflict_exists:
        cat = 'Internal conflict'
    else:
        cat, subcat = inheritance(
            *[example_files[family_rel] for family_rel in family_rels]
        )
    if args.create_qc_plots:
        for v in vcf_file_objs:
            if cat != '[Indel]' and v.has_curr_variant != 'absent':
                update_qc_counts(
                    qc_counts, v, cat, len([
                        v for v in vcf_file_objs if
                        v.has_curr_variant != 'absent'
                    ])
                )

    # Output each (REF, ALT, GT) combination at given position on chromosome
    if high_qual_present:
        for x in files_per_variant:
            output_ref_alt_gt(
                vcf_file_objs, qc_counts, overall_counts, args, out_all, x,
                files_per_variant, cat, subcat
            )


def output_ref_alt_gt(
    vcf_file_objs, qc_counts, overall_counts, args, out_all, ref_alt_gt,
    files_per_variant, cat, subcat
):
    """Process a variant uniquely identified by combination of
    CHROM, POS, REF, ALT, and GT: Write to file, and update
    overall_counts.
    """
    # List of VCF files that have given REF, ALT, GT
    ref_alt_gt_files = files_per_variant[ref_alt_gt]
    chrom = ref_alt_gt_files[0].curr_rec.CHROM
    pos = ref_alt_gt_files[0].curr_rec.POS
    ref, alt, gt = ref_alt_gt

    qc_averages = average_qc_values(ref_alt_gt_files)
    percent_present = 100.0 * len(ref_alt_gt_files) / len(vcf_file_objs)

    # Output variant to file
    variant_in_file_str = ''
    for v in vcf_file_objs:
        if (
            v.has_curr_variant != 'absent' and ref == str(v.curr_rec.REF) and
            alt == str(v.curr_rec.ALT[0]) and gt == v.curr_call['GT']
        ):
            variant_in_file_str += v.has_curr_variant
        else:
            variant_in_file_str += 'absent'
        variant_in_file_str += '\t'
    variant_in_file_str = variant_in_file_str[:-1]

    out_all.write(
        '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f'
        '\t%.2f\n' %
        (
            chrom, str(pos), ref, alt, gt, cat, subcat,
            variant_in_file_str, percent_present, qc_averages['qd'],
            qc_averages['dp'], qc_averages['qual'], qc_averages['ad0'],
            qc_averages['ad1']
        )
    )

    # Update overall counts
    overall_counts['Pass QC'] += 1
    if True in [
        v.curr_rec.is_indel for v in vcf_file_objs if
        v.has_curr_variant != 'absent'
    ]:
        overall_counts['Indel'] += 1
    else:
        overall_counts['No indel'] += 1
        overall_counts[cat] += 1
        if cat == 'Inherited':
            overall_counts[subcat] += 1


def internal_conflict_exists(vcf_file_objs, family_rel, files_per_variant):
    """Returns whether there is an internal conflict within files with
    given family_rel."""
    conflict_exists = False
    first_instance = None
    for x in files_per_variant:  # x == (REF, ALT, GT)
        for v in files_per_variant[x]:
            if v.family_rel == family_rel:
                if v.has_curr_variant == 'high_qual':
                    if first_instance in [None, x]:
                        first_instance = x
                    else:
                        conflict_exists = True
    curr_variants = [
        v.has_curr_variant for v in vcf_file_objs if
        v.family_rel == family_rel
    ]
    if 'absent' in curr_variants and 'high_qual' in curr_variants:
        conflict_exists = True
    return conflict_exists


def update_qc_counts(qc_counts, vcf_file_obj, cat, scale_factor):
    """Update qc_counts dict."""
    record = vcf_file_obj.curr_rec
    call = vcf_file_obj.curr_call
    qc_vals = {'dp': call['DP'], 'qd': record.INFO['QD'], 'qual': record.QUAL}
    for qc_cat in qc_vals:
        qc_val = qc_vals[qc_cat]
        intervals = qc_counts[qc_cat]['All']
        for x in intervals:
            if qc_val >= x:
                qc_counts[qc_cat][cat][x] += 1.0 / scale_factor


def average_qc_values(vcf_file_objs):
    """Returns dict containing average value for each QC category."""
    qc_avgs = {}
    qc_cats = [z for z in vcf_file_objs[0].curr_qc_values]
    for cat in qc_cats:
        qc_avgs[cat] = []
    for v in vcf_file_objs:
        if v.has_curr_variant == 'high_qual':
            for cat in v.curr_qc_values:
                qc_avgs[cat].append(v.curr_qc_values[cat])

    for cat in qc_cats:
        qc_avgs[cat] = (
            (float(sum(qc_avgs[cat])) / len(qc_avgs[cat])) if qc_avgs[cat]
            else -1
        )
    return qc_avgs


def min_chromosome(chroms):
    """Args: Chromosome names ['chrN', 'chrP', ...]. Returns name of
    chromosome which would appear first in ordered file.
    """
    alpha_chrom_order = ['X', 'Y', 'U', 'M']
    chroms_processed = list(chroms)
    for i in xrange(len(chroms_processed)):
        chroms_processed[i] = (
            chroms_processed[i][chroms_processed[i].index('chr')+3:]
        )
        try:
            chroms_processed[i] = int(chroms_processed[i])
        except ValueError:
            pass

    if any([isinstance(z, int) for z in chroms_processed]):
        return 'chr%s' % min(
            [y for y in chroms_processed if isinstance(y, int)]
        )
    else:
        first_letter = alpha_chrom_order[
            min([alpha_chrom_order.index(z[0]) for z in chroms_processed])
        ]
        return 'chr%s' % min(
            [y for y in chroms_processed if y[0] == first_letter]
        )


def print_dict_list(
    out_f, table_title, dict_list, cats, header_list, table_message
):
    """Create HTML file and print a list of category:value
    dictionaries as an HTML table.
    """
    if len(dict_list) == 0:
        return
    out_f.write('<h2>%s</h2>' % table_title)
    out_f.write(table_message)
    out_f.write(
        '<table border="1" style="table-layout:fixed; width:'
        '%dpx; word-wrap:break-word"><tr><th></th>' %
        (175 * (len(dict_list) + 1))
    )
    for x in header_list:
        out_f.write('<th>%s</th>' % x)
    out_f.write('</tr>')
    for cat in cats:
        out_f.write('<tr><th>%s</th>' % cat)
        for d in dict_list:
            out_f.write('<td>%s</td>' % (d[cat]))
        out_f.write('</tr>')
    out_f.write('</table>')


def normalize_counts(categories, counts):
    """Add normalized values to counts dict. Used for displaying
    HTML table.
    """
    if len(counts) == 0:
        return []
    min_index = 0
    # Normalized table entries at least flag_dist
    # distance from 1 are flagged (colored).
    flag_dist = .05
    counts_norm = [dict(z) for z in counts]
    for cat in categories:
        min_val = counts[min_index][cat]
        for f in xrange(len(counts)):
            if all(
                z[cat] == 0 and (not isinstance(z[cat], bool))
                for z in counts
            ):
                counts_norm[f][cat] = str(counts[f][cat]) + ' (%.2f)' % 1
            else:
                if (
                    isinstance(counts[f][cat], int) or
                    isinstance(counts[f][cat], float) and not
                    isinstance(counts[f][cat], bool) and
                    min_val != 0
                ):
                    counts_norm[f][cat] = (
                        counts[f][cat] / float(min_val) *
                        (
                            counts[min_index]['High-quality variants'] /
                            float(counts[f]['High-quality variants'])
                        )
                    )

                    if abs(1.0 - counts_norm[f][cat]) < flag_dist:
                        if isinstance(counts[f][cat], float):
                            counts_norm[f][cat] = '%.2f (%.2f)' % (
                                counts[f][cat], counts_norm[f][cat]
                            )
                        else:
                            counts_norm[f][cat] = '%d (%.2f)' % (
                                counts[f][cat], counts_norm[f][cat]
                            )
                    else:
                        if isinstance(counts[f][cat], float):
                            counts_norm[f][cat] = (
                                '%.2f (<font color = "blue">%.2f</font>)' % (
                                    counts[f][cat], counts_norm[f][cat]
                                )
                            )
                        else:
                            counts_norm[f][cat] = (
                                '%d (<font color = "blue">%.2f</font>)' % (
                                    counts[f][cat], counts_norm[f][cat]
                                )
                            )
    return counts_norm


def output_qc_plot(
    qc_categories, qc_param_label, qc_counts, out_filename, intervals
):
    """Create QC plot using pyplot. Save plot to out_filename."""
    plt.clf()
    for categ in qc_categories:
        if categ != 'All' and qc_counts[categ][0]:
            plt.plot(intervals, [
                (
                    1.0 - (float(qc_counts[categ][n]) / qc_counts[categ][0])
                ) for n in intervals
            ], label=categ)
    plt.axis([0, max(intervals), 0, 1.05])
    plt.xlabel(qc_param_label)
    plt.ylabel('Fraction removed')
    plt.title('%s as a QC parameter' % qc_param_label)
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * .75, box.height])
    ax.legend(
        [z for z in qc_categories if z != 'All'], loc='center left',
        bbox_to_anchor=(1, .9), prop={'size': 11.5}
    )
    plt.savefig(out_filename)


def inheritance(parent1_file, parent2_file, child_file):
    """Child's relationship to parents. Args: VCF file objects.
    Returns inheritance category and subcategory.
    VCF file obj representative of category. Assumes no internal
    conflicts
    """
    de_novo_label = 'De novo'
    invalid_label = 'Invalid'
    inherited_label = 'Inherited'

    cat = inherited_label
    subcat = 'NA'

    d = {
        'parent1': {'file': parent1_file}, 'parent2': {'file': parent2_file},
        'child': {'file': child_file}
    }
    for family_rel in d:
        if d[family_rel]['file']:
            d[family_rel]['ref'] = str(d[family_rel]['file'].curr_rec.REF)
            d[family_rel]['alt'] = str(d[family_rel]['file'].curr_rec.ALT[0])
            if d[family_rel]['file'].curr_call.gt_type == 1:
                d[family_rel]['gt'] = 'HET'
            elif d[family_rel]['file'].curr_call.gt_type == 2:
                d[family_rel]['gt'] = 'HOM-ALT'
            else:
                d[family_rel]['gt'] = None
        else:
            d[family_rel]['ref'] = None
            d[family_rel]['alt'] = None
            d[family_rel]['gt'] = 'HOM-REF'

    for x in d:
        if d[x]['gt'] != 'HOM-REF':
            chrom = d[x]['file'].curr_rec.CHROM
            pos = d[x]['file'].curr_rec.POS
    refs, alts = [], []
    for family_rel in d:
        v = d[family_rel]['file']
        if (v is not None and v.has_curr_variant == 'high_qual'):
            alts.append(str(v.curr_rec.ALT[0]))
            refs.append(str(v.curr_rec.REF))

    # Process variants for which there is a REF/ ALT/ GT conflict
    # between Parent 1 and Parent 2, or Parent 1 and Child, or
    # Parent 2 and Child
    if len(set(refs)) > 1:
        cat = '?'
        subcat = 'Child/ parent REFs don\'t match'
    elif len(set(alts)) > 1:
        if d['child']['gt'] == 'HOM-REF':
            if (
                d['parent1']['gt'] in ['HOM-REF', 'HET'] and
                d['parent2']['gt'] in ['HOM-REF', 'HET']
            ):
                cat = inherited_label
                subcat = 'Parents %s, %s; child HOM-REF' % ordered_gts(
                    d['parent1']['gt'], d['parent2']['gt']
                )
            else:
                cat = invalid_label
                subcat = 'Child HOM-REF; parents not both HOM-REF/ HET'
        elif d['child']['gt'] == 'HET':
            if (
                (
                    d['child']['alt'] == d['parent1']['alt'] and
                    d['parent1']['gt'] in ['HET', 'HOM-ALT'] and
                    d['parent2']['gt'] in ['HOM-REF', 'HET']
                ) or
                (
                    d['child']['alt'] == d['parent2']['alt'] and
                    d['parent2']['gt'] in ['HET', 'HOM-ALT'] and
                    d['parent1']['gt'] in ['HOM-REF', 'HET']
                )
            ):
                cat = inherited_label
                subcat = 'Parents %s, %s; child HET' % ordered_gts(
                    d['parent1']['gt'], d['parent2']['gt']
                )
            elif (
                d['parent1']['gt'] == 'HOM-REF' and
                d['parent2']['gt'] == 'HOM-REF'
            ):
                cat = de_novo_label
            else:
                cat = invalid_label
                subcat = (
                    'Child HET; not the case that one parent has '
                    'child ALT and other parent has REF'
                )
        elif d['child']['gt'] == 'HOM-ALT':
            if (
                d['parent1']['gt'] == 'HOM-REF' and
                d['parent2']['gt'] == 'HOM-REF'
            ):
                cat = de_novo_label
            else:
                cat = invalid_label
                subcat = (
                    'Child HOM-ALT; parents not both HET/ HOM-ALT for '
                    'child ALT allele'
                )
        else:
            cat = invalid_label
            subcat = 'Unknown/ invalid child GT'

    else:
        if d['child']['gt'] == 'HOM-REF':
            if (
                d['parent1']['gt'] not in ['HOM-REF', 'HET'] or
                d['parent2']['gt'] not in ['HOM-REF', 'HET']
            ):
                cat = invalid_label
                subcat = 'Child HOM-REF; parents not both HOM-REF/ HET'
        elif d['child']['gt'] == 'HET':
            if (
                d['parent1']['gt'] == 'HOM-ALT' and
                d['parent2']['gt'] == 'HOM-ALT'
            ):
                cat = invalid_label
                subcat = 'Child HET; parents both HOM-ALT'
        elif d['child']['gt'] == 'HOM-ALT':
            if not (
                d['parent1']['gt'] in ['HET', 'HOM-ALT'] and
                d['parent2']['gt'] in ['HET', 'HOM-ALT']
            ):
                cat = invalid_label
                subcat = 'Child HOM-ALT; parents not both HET/ HOM-ALT'
        else:
            cat = invalid_label
            subcat = 'Unknown/ invalid child GT'
        if cat != inherited_label:
            return (cat, subcat)

        if (
            d['parent1']['gt'] == 'HOM-REF' and
            d['parent2']['gt'] == 'HOM-REF'
        ):
            if d['child']['gt'] == 'HOM-REF':
                subcat = 'Parents HOM-REF, HOM-REF; child HOM-REF'
            else:
                cat = de_novo_label
        elif (
            (
                d['parent1']['gt'] == 'HOM-REF' and
                d['parent2']['gt'] == 'HET'
            ) or (
                d['parent1']['gt'] == 'HET' and
                d['parent2']['gt'] == 'HOM-REF'
            )
        ):
            if d['child']['gt'] == 'HOM-REF':
                subcat = 'Parents HOM-REF, HET; child HOM-REF'
            elif d['child']['gt'] == 'HET':
                subcat = 'Parents HOM-REF, HET; child HET'
            else:
                cat = invalid_label
                subcat = 'Parents HOM-REF, HET; child not HOM-REF/ HET'
        elif (
            (
                d['parent1']['gt'] == 'HOM-REF' and
                d['parent2']['gt'] == 'HOM-ALT'
            ) or (
                d['parent1']['gt'] == 'HOM-ALT' and
                d['parent2']['gt'] == 'HOM-REF')
        ):
            if d['child']['gt'] == 'HET':
                subcat = 'Parents HOM-REF, HOM-ALT; child HET'
            else:
                cat = invalid_label
                subcat = 'Parents HOM-REF, HOM-ALT; child not HET'
        elif (
            d['parent1']['gt'] == 'HET' and
            d['parent2']['gt'] == 'HET'
        ):
            if d['child']['gt'] == 'HOM-REF':
                subcat = 'Parents HET, HET; child HOM-REF'
            elif d['child']['gt'] == 'HET':
                subcat = 'Parents HET, HET; child HET'
            elif d['child']['gt'] == 'HOM-ALT':
                subcat = 'Parents HET, HET; child HOM-ALT'
            else:
                cat = invalid_label
                subcat = 'Unknown/ invalid child GT'
        elif (
            (
                d['parent1']['gt'] == 'HET' and
                d['parent2']['gt'] == 'HOM-ALT'
            ) or (
                d['parent1']['gt'] == 'HOM-ALT' and
                d['parent2']['gt'] == 'HET'
            )
        ):
            if d['child']['gt'] == 'HET':
                subcat = 'Parents HET, HOM-ALT; child HET'
            elif d['child']['gt'] == 'HOM-ALT':
                subcat = 'Parents HET, HOM-ALT; child HOM-ALT'
            else:
                cat = invalid_label
                subcat = 'Parents HET, HOM-ALT; child not HET/ HOM-ALT'
        elif (
            d['parent1']['gt'] == 'HOM-ALT' and
            d['parent2']['gt'] == 'HOM-ALT'
        ):
            if d['child']['gt'] == 'HOM-ALT':
                subcat = 'Parents HOM-ALT, HOM-ALT; child HOM-ALT'
            else:
                cat = invalid_label
                subcat = 'Parents HOM-ALT, HOM-ALT; child not HOM-ALT'
        else:
            cat = invalid_label
            subcat = 'Unknown/ invalid parent GT(s)'

    return (cat, subcat)


def ordered_gts(gt1, gt2):
    """Sorts two genotypes in ascending order"""
    gts = ['HOM-REF', 'HET', 'HOM-ALT']
    if gt1 not in gts or gt2 not in gts:
        return (None, None)
    return ((gt1, gt2) if (gts.index(gt1) < gts.index(gt2)) else (gt2, gt1))


def predict_sex(
    total_x_dp, total_y_dp, total_x_gt, total_y_gt, count_x, count_y
):
    """Predict whether individual is male or female based on
    X and Y chromosome variants. GT ranges from 0 (HET) to 1 (HOM-ALT).
    """
    avg_x_dp, avg_x_gt = None, None if (count_x == 0) else (
        float(total_x_dp) / count_x, float(total_x_gt) / count_x
    )
    avg_y_dp, avg_y_gt = None, None if (count_y == 0) else (
        float(total_y_dp) / count_y, float(total_y_gt) / count_y
    )

    female_score, male_score = 0, 0
    if avg_x_dp and avg_y_dp:
        if avg_x_dp > 2 * avg_y_dp:
            female_score += 1
        else:
            male_score += 1
    if avg_x_gt:
        if avg_x_gt < 0.65:
            female_score += 1
        else:
            male_score += 1
    if avg_y_gt:
        if avg_y_gt > 0.8:
            female_score += 1
        else:
            male_score += 1

    if female_score > male_score:
        return 'F'
    elif male_score > female_score:
        return 'M'
    else:
        return 'No prediction'

if __name__ == '__main__':
    sys.exit(main())
