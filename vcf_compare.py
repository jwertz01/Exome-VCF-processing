"""Compares VCF files. Produces statistics and variant comparison files."""

import os
import sys
import argparse
import matplotlib.pyplot as plt
import vcf


class VcfFile(object):
    """Contains info about VCF file, and current position in file."""
    def __init__(
        self, file_name=None, file_=None, counts=None, reader=None,
        curr_rec=None, next_rec=None, curr_call=None, eof=False,
        has_curr_variant='absent', curr_qc_values=None,
        is_rediscovery_file=False
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
        self.curr_qc_values = curr_qc_values  # DP, QUAL, etc. for current var
        self.is_rediscovery_file = is_rediscovery_file

    def __str__(self):
        return (
            'Filename: %s\nCounts: %s\nCurr rec: %s\nNext rec: %s\nEOF: %s'
            '\nHas curr variant: %s\nCurr QC values: %s\nIs rediscovery '
            'file: %s\n' %
            (
                self.file_name, self.counts, self.curr_rec, self.next_rec,
                self.eof, self.has_curr_variant, self.curr_qc_values,
                self.is_rediscovery_file
            )
        )


def main():
    parser = argparse.ArgumentParser()
    args = parse_arguments(parser)
    args.min_qc_values = {
        'dp': args.min_dp, 'qd': args.min_qd, 'qual': args.min_qual
    }
    validate_args(args, parser)

    # Per-file statistics, in the order in which they will be displayed
    # in HTML table. (Categories starting with a capital letter are
    # displayed.) Stored in counts dict in VcfFile object.
    counts_per_file_categs = [
        'Total variants', 'Average DP', 'Average QUAL', 'average QD',
        'total_dp', 'count_dp', 'total_qual', 'count_qual', 'total_qd',
        'count_qd', 'High-quality variants', 'Homozygous reference', 'SNP',
        'Heterozygous', 'Homozygous alternate', 'Indel', 'Missing',
        'a>c', 'a>g', 'a>t', 'c>a', 'c>g', 'c>t', 'g>a', 'g>c', 'g>t', 't>a',
        't>c', 't>g', 'Transitions/ transversions',
        'Variants absent from this file present in other files',
        'Variants present that conflict with variants present in other files',
        'Variants unique to this file'
    ]
    for f_name in args.rediscovery_files:
        counts_per_file_categs.append(
            'Rediscovery rate (non-indels): %s' % f_name
        )
        counts_per_file_categs.append('rediscovery_count_%s' % f_name)

    # Categories displayed in flowchart in HTML file
    overall_counts_categs = [
        'Pass QC', 'All agree', 'Conflicts', 'Conflicts absent',
        'Conflicts multiple', 'Conflicts multiple indel',
        'Conflicts multiple no indel', 'Low QC', 'Agree GT disagree ALT',
        'Disagree GT disagree ALT', 'Disagree GT agree ALT', 'Predict HET',
        'Predict_HOM_ALT', 'No prediction'
    ]
    overall_counts = {}

    # Categories displayed in QC plots
    qc_counts_categs = [
        'All', 'HET all agree', 'HET disagree', 'HOM all agree',
        'HOM disagree'
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
    empty_files = []
    for f_name in args.vcf_file_names + args.rediscovery_files:
        in_f = open(f_name)
        reader = vcf.Reader(in_f)
        f_base_name = os.path.basename(f_name)
        if len(reader.samples) > 1:
            multi_sample_files.append(f_base_name)
            in_f.close()
        else:
            v = VcfFile(
                file_name=f_base_name, file_=in_f, counts={}, reader=reader,
                curr_qc_values={}
            )
            try:
                v.next_rec = reader.next()
            except StopIteration:
                v.eof = True
            for cat in counts_per_file_categs:
                v.counts[cat] = 0
            if v.reader and v.next_rec:
                vcf_file_objs.append(v)
            else:
                empty_files.append(f_base_name)
            if f_name in args.rediscovery_files:
                v.is_rediscovery_file = True
    if args.create_qc_plots:
        for i in qc_intervals:
            qc_counts[i] = {}
            for j in qc_counts_categs:
                qc_counts[i][j] = {}
                for k in qc_intervals[i]:
                    qc_counts[i][j][k] = 0
    for x in overall_counts_categs:
        overall_counts[x] = 0

    # Iterate over vcf files, write variant comparison output files,
    # fill in counts
    header = (
        'CHROM\tPOS\tREF\tALT\tGT\t%s\t%%Present\tAvg_QD\tAvg_DP\t'
        'Avg_QUAL\tAvg_AD0\tAvg_AD1' %
        ('\t'.join([v.file_name for v in vcf_file_objs]))
    )
    out_files = []
    with \
            open(args.out_all_variants, 'w') as out_all, \
            open(args.out_absent, 'w') as out_absent, \
            open(args.out_multiple, 'w') as out_multiple:
        out_all.write(header + '\n')
        out_absent.write(header + '\n')
        out_multiple.write(header + '\tSummary\tConclusion\n')
        out_files = {
            'all': out_all, 'absent': out_absent, 'multiple': out_multiple
        }

        compare_variants(
            vcf_file_objs, qc_counts, overall_counts, args, out_files
        )

    for v in vcf_file_objs:
        v.file_.close()
        if v.is_rediscovery_file:
            continue

        v.counts['Average DP'] = divide(
            v.counts['total_dp'], v.counts['count_dp'], -1.0
        )
        v.counts['Average QUAL'] = divide(
            v.counts['total_qual'], v.counts['count_qual'], -1.0
        )
        v.counts['Average QD'] = divide(
            v.counts['total_qd'], v.counts['count_qd'], -1.0
        )
        for y in [z for z in vcf_file_objs if z.is_rediscovery_file]:
            v.counts[
                'Rediscovery rate (non-indels): %s' % y.file_name
            ] = divide(
                v.counts['rediscovery_count_%s' % y.file_name],
                v.counts['SNP'], -1
            )
        v.counts['Transitions/ transversions'] = ts_tv_ratio(v.counts)

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
        counts_per_file_categs,
        [v.counts for v in vcf_file_objs if not v.is_rediscovery_file]
    )
    with open(args.out_stats, 'w') as out_f:
        page_title = 'VCF Comparison'
        out_f.write(
            '<!DOCTYPE html><html><head><title>%s</title></head><body>'
            '<h1>%s</h1><br />' % (page_title, page_title)
        )

        categs_pre_filter = ['Total variants', 'Average DP', 'Average QUAL']

        # per unique combination of CHROM and POS
        categs_post_filter_1 = [
            'High-quality variants', 'Homozygous reference', 'SNP',
            'Heterozygous', 'Homozygous alternate', 'Indel', 'Missing',
            'Transitions/ transversions'
        ]
        for f_name in args.rediscovery_files:
            categs_post_filter_1.append(
                'Rediscovery rate (non-indels): %s' % f_name
            )

        # per unique combination of CHROM, POS, REF, ALT, GT
        categs_post_filter_2 = [
            'Variants absent from this file present in other files',
            'Variants present that conflict with variants present '
            'in other files',
            'Variants unique to this file'
        ]

        print_dict_list(
            out_f, 'Statistics per file, all variants,<br />per unique '
            'combination of CHROM and POS:', counts_norm, categs_pre_filter,
            [v.file_name for v in vcf_file_objs if not v.is_rediscovery_file]
        )
        print_dict_list(
            out_f, 'Statistics per file, high-quality variants only,<br />'
            'per unique combination of CHROM and POS:', counts_norm,
            categs_post_filter_1,
            [v.file_name for v in vcf_file_objs if not v.is_rediscovery_file]
        )
        print_dict_list(
            out_f, 'Statistics per file, high-quality variants only,<br />'
            'per unique combination of CHROM, POS, REF, ALT, and GT:',
            counts_norm, categs_post_filter_2,
            [v.file_name for v in vcf_file_objs if not v.is_rediscovery_file]
        )
        if len(multi_sample_files) > 0:
            out_f.write(
                '<p><b>Multi-sample files excluded from ' +
                'analysis: </b>%s</p>' % (', '.join(multi_sample_files)))
        if len(empty_files) > 0:
            out_f.write(
                '<p><b>Empty/ invalid files excluded from ' +
                'analysis: </b>%s</p>' % (', '.join(empty_files)))
        out_f.write(
            '<h2>Overall statistics, high-quality variants only,<br />'
            'per unique combination of CHROM, POS, REF, ALT, and GT:</h2>'
        )
        out_f.write(
            '<p>Each category (except Pass QC) is a subcategory of the '
            'rightmost category directly above it.</p>'
        )
        out_f.write('<table><tr>')
        for x in [
            z for z in overall_counts_categs if 'predict' not in z.lower()
        ]:
            if x in [
                'All agree', 'Conflicts absent',
                'Conflicts multiple indel', 'Low QC'
            ]:
                out_f.write('</tr><tr>')
            out_f.write(
                '<td style=padding-right:2em><b>%s</b>: %d</td>' %
                (x, overall_counts[x])
            )
        out_f.write('</tr></table>')
        out_f.write('</body></html>')


def divide(num, denom, default):
    """Return numerator divided by denominator"""
    if denom > 0:
        return float(num) / denom
    else:
        return default


def parse_arguments(parser):
    """Add arguments to parser. Return parsed arguments."""
    parser.add_argument(
        'vcf_file_names', nargs='+',
        help='Names of VCF files to be processed (required)'
    )
    parser.add_argument(
        '--create_qc_plots', action='store_true',
        help='Whether to output QC plots. Takes no arguments. '
             '(Default: Does not output plots.)'
    )
    parser.add_argument(
        '--rediscovery_files', nargs='+', default=[],
        help='Names of VCF files to be used for '
             'calculating rediscovery rate. (Default: Does not '
             'calculate rediscovery rate.)'
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
        '--out_absent', default='conflicts_absent.txt',
        help='Name of output file containing variants such that at least one '
             'file has high-quality variant and another file does not have '
             'variant at CHROM and POS (Default: conflicts_absent.txt)'
    )
    parser.add_argument(
        '--out_multiple',
        default='conflicts_multiple.txt',
        help='Name of output file containing variants such that at least one '
             'file has high-quality variant and another file has different '
             'REF/ALT/GT at CHROM and POS. Excludes indels. (Default: '
             'conflicts_multiple.txt)'
    )
    parser.add_argument(
        '--min_qd', type=int, default=0,
        help='Minimum quality by depth for inclusion in variant '
             'output files. (Default: 0)'
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


def validate_args(args, parser):
    """Validate command line arguments."""
    message = ''
    for x in args.vcf_file_names + args.rediscovery_files:
        message += check_extension(x, '.vcf')
    message += check_extension(args.out_stats, '.html')
    message += check_extension(args.out_all_variants, '.txt')
    message += check_extension(args.out_absent, '.txt')
    message += check_extension(args.out_multiple, '.txt')

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
    vcf_file_objs, qc_counts, overall_counts, args, out_files
):
    """Iterate through all VCF files simultaneously. Write per-file
    variant info to output files and update overall statistics
    (counts).
    """
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
        min_pos = min([
            v.next_rec.POS for v in vcf_file_objs if (
                parse_chromosome(v.next_rec.CHROM) == parse_chromosome(
                    min_chrom
                ) and (not v.eof)
            )
        ])
        for v in vcf_file_objs:
            if (v.reader is None) or (v.next_rec is None) or (v.counts == {}):
                continue
            v.curr_rec = v.next_rec
            v.curr_call = None
            record = v.next_rec
            if v.reader.samples:
                call = record.genotype(v.reader.samples[0])
            else:
                call = None
            if (
                parse_chromosome(record.CHROM) == parse_chromosome(min_chrom)
                and record.POS == min_pos
                and not v.eof
            ):
                if call:
                    assign_qc_values(v.curr_qc_values, record, call)
                    update_counts_all(v.counts, v.curr_qc_values)

                    # count as high-quality if qc values are unknown
                    is_high_qual = True
                    for qc_param in args.min_qc_values:
                        if (
                            v.curr_qc_values[qc_param] >= 0 and (
                                v.curr_qc_values[qc_param] <
                                args.min_qc_values[qc_param]
                            )
                        ):
                            is_high_qual = False
                    v.curr_call = call

                    if is_high_qual:
                        v.has_curr_variant = 'high_qual'
                        update_counts_high_qual(v.counts, record, call)

                        for y in vcf_file_objs:
                            if (
                                y.is_rediscovery_file and
                                parse_chromosome(
                                    y.next_rec.CHROM
                                ) == parse_chromosome(min_chrom) and
                                y.next_rec.POS == min_pos and
                                not record.is_indel
                            ):
                                v.counts[
                                    'rediscovery_count_%s' % y.file_name
                                ] += 1
                    else:
                        v.has_curr_variant = 'low_qual'
                else:
                    v.has_curr_variant = 'present'
                    v.counts['Total variants'] += 1
                    v.counts['High-quality variants'] += 1
                try:
                    v.next_rec = v.reader.next()

                except StopIteration:
                    v.eof = True
            else:
                v.has_curr_variant = 'absent'

        # Write single variant (at min_chrom and min_pos) to output files,
        # update counts
        process_variant(
            vcf_file_objs, qc_counts, overall_counts, args, out_files
        )


def assign_qc_values(qc_values_dict, record, call):
    """Assign QC values to dict, given current record and call."""
    try:
        qc_values_dict['qd'] = record.INFO['QD']
    except (KeyError, AttributeError):
        qc_values_dict['qd'] = -1.0
    try:
        qc_values_dict['dp'] = call['DP']
    except (KeyError, AttributeError):
        try:
            qc_values_dict['dp'] = record.INFO['DP']
        except (KeyError, AttributeError):
            qc_values_dict['dp'] = -1.0
    try:
        qc_values_dict['qual'] = record.QUAL
    except (KeyError, AttributeError):
        qc_values_dict['qual'] = -1.0
    try:
        qc_values_dict['ad0'] = call['AD'][0]
    except (KeyError, AttributeError):
        qc_values_dict['ad0'] = -1
    try:
        qc_values_dict['ad1'] = call['AD'][1]
    except (KeyError, AttributeError):
        qc_values_dict['ad1'] = -1


def update_counts_all(counts, qc_values):
    """Update counts dict for a variant."""
    counts['Total variants'] += 1
    if qc_values['dp'] >= 0:
        counts['total_dp'] += qc_values['dp']
        counts['count_dp'] += 1
    if qc_values['qual'] >= 0:
        max_qual = 1000000  # Exclude outliers
        if qc_values['qual'] > max_qual:
            counts['total_qual'] += max_qual
        else:
            counts['total_qual'] += qc_values['qual']
        counts['count_qual'] += 1
    if qc_values['qd'] >= 0:
        counts['total_qd'] += qc_values['qd']
        counts['count_qd'] += 1


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
        snp_type = '%s>%s' % (
            record.alleles[0],
            record.alleles[int(max(call.gt_alleles))]
        )
        if snp_type.lower() in counts:
            counts[snp_type.lower()] += 1
        else:
            counts[snp_type.lower()] = 1


def process_variant(
    vcf_file_objs, qc_counts, overall_counts, args, out_files
):
    """Process a variant uniquely identified by CHROM and POS.
    Identify and process sub-variants."""
    high_qual_present = any([
        v.has_curr_variant == 'high_qual' for v in vcf_file_objs
    ])
    files_per_variant = {}     # (REF, ALT, GT): list of VCF file objs

    if high_qual_present:
        for v in vcf_file_objs:
            if v.has_curr_variant in ['low_qual', 'high_qual']:
                chrom = parse_chromosome(v.curr_rec.CHROM)
                pos = v.curr_rec.POS
        for v in vcf_file_objs:
            if v.has_curr_variant in ['low_qual', 'high_qual']:
                call = v.curr_call
                record = v.curr_rec
                ref = str(record.REF)
                alt = str(record.ALT[0])
                gt = call['GT']

                if (ref, alt, gt) in files_per_variant:
                    files_per_variant[(ref, alt, gt)].append(v)
                else:
                    files_per_variant[(ref, alt, gt)] = [v]

        # Output each (REF, ALT, GT) combination at given position
        # on chromosome
        for x in files_per_variant:
            output_ref_alt_gt(
                x, out_files, vcf_file_objs, overall_counts, qc_counts,
                args, files_per_variant
            )

    # Update plot data
    if args.create_qc_plots:
        multiple_variants = [
            'high_qual' in [v.has_curr_variant for v in files_per_variant[x]]
            for x in files_per_variant
        ].count(True) > 1

        all_agree = (
            all(
                v.has_curr_variant != 'absent' for v in
                vcf_file_objs
            ) and (not multiple_variants)
        )

        for v in vcf_file_objs:
            if v.has_curr_variant not in ['low_qual', 'high_qual']:
                continue
            gt_type = v.curr_call.gt_type
            cat = 'undefined'
            if all_agree:
                if gt_type == 1:
                    cat = 'HET all agree'
                elif gt_type == 2:
                    cat = 'HOM all agree'
            else:
                if gt_type == 1:
                    cat = 'HET disagree'
                elif gt_type == 2:
                    cat = 'HOM disagree'
            if cat != 'undefined':
                update_qc_counts(
                    qc_counts, v, cat, len([
                        v for v in vcf_file_objs if
                        v.has_curr_variant in ['low_qual', 'high_qual']
                    ])
                )


def output_ref_alt_gt(
    ref_alt_gt, out_files, vcf_file_objs, overall_counts, qc_counts,
    args, files_per_variant
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
    percent_present = divide(
        100 * len(ref_alt_gt_files),
        len([z for z in vcf_file_objs if not z.is_rediscovery_file]), -1.0
    )

    multiple_variants = [
        'high_qual' in [
            v.has_curr_variant for v in files_per_variant[x] if not
            v.is_rediscovery_file
        ] for x in files_per_variant
    ].count(True) > 1

    variant_in_file_str = ''
    for v in vcf_file_objs:
        if (
                v.has_curr_variant != 'absent' and
                ref == str(v.curr_rec.REF) and
                alt == str(v.curr_rec.ALT[0]) and
                ((not v.curr_call) or (gt == v.curr_call['GT']))
        ):
            variant_in_file_str += v.has_curr_variant
        else:
            variant_in_file_str += 'absent'
        variant_in_file_str += '\t'
    variant_in_file_str = variant_in_file_str[:-1]

    for v in ref_alt_gt_files:
        if (
            v.has_curr_variant == 'high_qual' and percent_present != 100
            and multiple_variants
        ):
            v.counts[
                'Variants present that conflict with variants present '
                'in other files'
            ] += 1
        if (
            v.has_curr_variant == 'high_qual' and len(
                [z for z in ref_alt_gt_files if not z.is_rediscovery_file]
            ) == 1
        ):
            v.counts['Variants unique to this file'] += 1

    for v in vcf_file_objs:
        if v not in ref_alt_gt_files and any(
            [not y.is_rediscovery_file for y in ref_alt_gt_files]
        ):
            v.counts[
                'Variants absent from this file present in other files'
            ] += 1

    out_str = (
        '%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n'
        % (
            chrom, str(pos), ref, alt, gt, variant_in_file_str,
            percent_present, qc_averages['qd'], qc_averages['dp'],
            qc_averages['qual'], qc_averages['ad0'], qc_averages['ad1']
        )
    )

    # Output variant to files and update overall counts.
    if any([(not z.is_rediscovery_file) for z in ref_alt_gt_files]):
        overall_counts['Pass QC'] += 1
        out_files['all'].write(out_str)
        if percent_present == 100:
            overall_counts['All agree'] += 1
        else:
            overall_counts['Conflicts'] += 1
            if multiple_variants:
                overall_counts['Conflicts multiple'] += 1
                if True in [
                    v.curr_rec.is_indel for v in vcf_file_objs if
                    v.has_curr_variant != 'absent'
                ]:
                    overall_counts['Conflicts multiple indel'] += 1
                else:
                    overall_counts['Conflicts multiple no indel'] += 1
                    out_files['multiple'].write(out_str[:-1])
                    write_summary_conclusion(
                        files_per_variant, out_files['multiple'],
                        overall_counts
                    )
            else:
                overall_counts['Conflicts absent'] += 1
                out_files['absent'].write(out_str)


def update_qc_counts(qc_counts, vcf_file_obj, cat, scale_factor):
    """Update qc_counts dict."""
    record = vcf_file_obj.curr_rec
    call = vcf_file_obj.curr_call
    for qc_cat in qc_counts:
        qc_value = vcf_file_obj.curr_qc_values[qc_cat]
        intervals = qc_counts[qc_cat]['All']
        for x in intervals:
            if isinstance(scale_factor, int) and scale_factor > 0:
                if qc_value >= x:
                    qc_counts[qc_cat][cat][x] += 1.0 / scale_factor
            else:
                raise ValueError('Scale factor must be positive integer.')


def average_qc_values(vcf_file_objs):
    """Return dict containing average value for each QC category."""
    qc_avgs = {}
    qc_cats = [z for z in vcf_file_objs[0].curr_qc_values]
    for cat in qc_cats:
        qc_avgs[cat] = []
    for v in vcf_file_objs:
        if v.has_curr_variant == 'high_qual':
            for cat in v.curr_qc_values:
                qc_avgs[cat].append(v.curr_qc_values[cat])

    for cat in qc_cats:
        if None in qc_avgs[cat]:
            qc_avgs[cat] = -1
        else:
            qc_avgs[cat] = divide(sum(qc_avgs[cat]), len(qc_avgs[cat]), -1)
    return qc_avgs


def write_summary_conclusion(files_per_variant, out_multiple, overall_counts):
    """Write variant info (summary and conclusion) to out_multiple file,
    and update overall_counts dict.
    """
    conflicts = conflict_categs(files_per_variant)
    common_gt = write_most_frequent(
        conflicts, files_per_variant, out_multiple
    )
    ad_ratio_categ = write_ad_ratio(
        files_per_variant, out_multiple, .1, .8, 10
    ) if ('GT' in conflicts) else 'other'
    max_qd = write_qd_summary(files_per_variant, out_multiple, 5, 10)
    write_conclusion(
        overall_counts, conflicts, ad_ratio_categ, common_gt, out_multiple,
        max_qd, 10
    )


def conflict_categs(files_per_variant):
    """Find whether REF/ ALT/ GT conflict is present."""
    conflicts = []
    if len(set([r for (r, a, g) in files_per_variant])) > 1:
        conflicts.append('REF')
    if len(set([a for (r, a, g) in files_per_variant])) > 1:
        conflicts.append('ALT')
    if len(set([g for (r, a, g) in files_per_variant])) > 1:
        conflicts.append('GT')
    return conflicts


def write_most_frequent(conflicts, files_per_variant, out_f):
    """Write REF/ ALT/ GT value that occurs most frequently in files,
    if conflict of that type is present.
    """
    plur = []
    common_gt = None
    if 'REF' in conflicts:
        common_ref = plurality(files_per_variant, 0)
        plur += common_ref
    if 'ALT' in conflicts:
        common_alt = plurality(files_per_variant, 1)
        plur += common_alt
    if 'GT' in conflicts:
        common_gt = plurality(files_per_variant, 2)
        plur += common_gt
    out_f.write(
        '\t%s conflict.' % ', '.join(conflicts) + (
            (' Most frequent: %s.' % ', '.join(plur)) if (plur != []) else ''
        )
    )
    return common_gt


def plurality(tuple_dict, n):
    """
        Args:
            tuple_dict (dict of (a, b, c, ...) : list of VcfFile objs)
            n (int): Position in tuple_dict key
        Returns:
            Element(s) of (a, b, c, ...)[n] that occurs most frequently
            in files
    """
    list_ = [
        [x[n]] * [
            y.has_curr_variant != 'absent' for y in tuple_dict[x]
        ].count(True)
        for x in tuple_dict
    ]
    list_ = [x for y in list_ for x in y]  # Merge sublists
    max_count = max([list_.count(x) for x in list_])
    plur = list(set([x for x in list_ if list_.count(x) == max_count]))
    return ([] if len(plur) == len(set(list_)) else plur)


def write_ad_ratio(
    files_per_variant, out_f, hom_alt_ratio, het_ratio, hom_ref_ratio
):
    """Predict whether variant is HOM-ALT/ HET/ HOM-REF, based on
    ad0/ad1 ratio.
    """
    total_ad0, total_ad1 = 0, 0
    for x in files_per_variant:
        for v in files_per_variant[x]:
            if v.has_curr_variant == 'high_qual':
                ad0 = v.curr_qc_values['ad0']
                ad1 = v.curr_qc_values['ad1']
                if ad0 >= 0 and ad1 >= 0:
                    total_ad0 += ad0
                    total_ad1 += ad1

    if total_ad0 >= 0:
        ad_ratio = divide(total_ad0, total_ad1, -1.0)
    else:
        ad_ratio = -1.0
    ad_ratio_categ = 'other'
    if ad_ratio < 0:
        pass
    elif ad_ratio <= hom_alt_ratio:
        ad_ratio_categ = 'HOM-ALT'
    elif ad_ratio > het_ratio and ad_ratio < divide(
        1, het_ratio, float('-inf')
    ):
        ad_ratio_categ = 'HET'
    elif (ad_ratio > hom_ref_ratio):
        ad_ratio_categ = 'HOM-REF'

    if ad_ratio_categ != 'other':
        out_f.write(' AD0/AD1 ratio indicates %s.' % ad_ratio_categ)
    return ad_ratio_categ


def write_qd_summary(files_per_variant, out_f, low_qd, borderline_qd):
    """Write whether QD is low/ borderline."""
    try:
        max_qd = max([
            average_qc_values(files_per_variant[x])['qd']
            for x in files_per_variant
        ])
    except KeyError:
        max_qd = float('inf')

    if max_qd < low_qd:
        out_f.write(' Low QD.')
    elif max_qd < borderline_qd:
        out_f.write(' Borderline QD.')
    return max_qd


def write_conclusion(
    overall_counts, conflicts, ad_ratio_categ, common_gt, out_f, max_qd,
    borderline_qd
):
    """Write QC summary, type(s) of conflict if present, HET/ HOM-ALT
    prediction. Update overall_counts."""
    if (max_qd < borderline_qd):
        overall_counts['Low QC'] += 1
        out_f.write('\tNeeds validation (low/borderline QC).')
    elif 'REF' in conflicts:
        out_f.write('\tREF conflict(?)')
    elif 'GT' in conflicts and ad_ratio_categ == 'HOM-REF':
        out_f.write('\tHOM-REF AD ratio(?)')
    elif 'GT' in conflicts and 'ALT' in conflicts:
        overall_counts['Disagree GT disagree ALT'] += 1
        out_f.write('\tDisagree GT, ALT.')
    elif 'GT' in conflicts:
        overall_counts['Disagree GT agree ALT'] += 1
        out_f.write('\tDisagree GT, agree ALT.')
    elif 'ALT' in conflicts:
        overall_counts['Agree GT disagree ALT'] += 1
        out_f.write('\tAgree GT, disagree ALT.')
    else:
        out_f.write('\tNo conflicts(?)')
    out_f.write('\n')


def min_chromosome(chroms):
    """Args: Chromosome names ['chrN', 'chrP', ...]. Return name of
    chromosome which would appear first in ordered file.
    """
    alpha_chrom_order = ['X', 'Y', 'U', 'M']
    chroms_list = list(chroms)
    for i in xrange(len(chroms_list)):
        chroms_list[i] = parse_chromosome(chroms_list[i])
    if any([isinstance(z, int) for z in chroms_list]):
        return 'chr%s' % min(
            [y for y in chroms_list if isinstance(y, int)]
        )
    else:
        first_letter = alpha_chrom_order[
            min([alpha_chrom_order.index(z[0]) for z in chroms_list])
        ]
        return 'chr%s' % min(
            [y for y in chroms_list if y[0] == first_letter]
        )


def parse_chromosome(chrom):
    """Parse chromosome name."""
    if chrom.startswith('chr'):
        chrom = chrom[3:]
    try:
        chrom = int(chrom)
    except ValueError:
        pass
    return chrom


def print_dict_list(out_f, table_title, dict_list, cats, header_list):
    """Create HTML file and print a list of category:value dictionaries
    as an HTML table.
    """
    if len(dict_list) == 0:
        return
    out_f.write('<h2>%s</h2>' % table_title)
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
    out_f.write('<br />')


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
                    (
                        isinstance(counts[f][cat], int) or
                        isinstance(counts[f][cat], float) and not
                        isinstance(counts[f][cat], bool)
                    ) and (counts[f]['High-quality variants'] > 0)
                    and (counts[f][cat] > 0) and (min_val > 0)
                ):
                    # If normalization not relative to number of variants
                    if (
                        cat in [
                            'Transitions/ transversions',
                            'High-quality variants'
                        ] or 'Rediscovery rate' in cat or 'Average' in cat
                    ):
                        counts_norm[f][cat] = counts[f][cat] / float(min_val)
                    else:
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
            plt.plot(
                intervals, [
                    (1.0 - (float(qc_counts[categ][n]) / qc_counts[categ][0]))
                    for n in intervals
                ], label=categ
            )
    plt.axis([0, max(intervals), 0, 1.05])
    plt.xlabel(qc_param_label)
    plt.ylabel('Fraction removed')
    plt.title('%s as a QC parameter' % qc_param_label)
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * .75, box.height])
    ax.legend(
        [z for z in qc_categories if z != 'All'], loc='center left',
        bbox_to_anchor=(1, .87))
    plt.savefig(out_filename)


def ts_tv_ratio(counts):
    """Calculate ratio of transitions / transversions."""
    ts = ['a>g', 'g>a', 'c>t', 't>c']
    tv = ['a>c', 'c>a', 'a>t', 't>a', 'c>g', 'g>c', 'g>t',  't>g']
    return divide(
        sum([counts[x] for x in ts]), sum([counts[y] for y in tv]), -1.0
    )


if __name__ == '__main__':
    sys.exit(main())
