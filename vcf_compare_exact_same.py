"""Compares VCF files that are expected to be exactly the same.
Produces statistics and variant comparison files."""

import os
import sys
import argparse
import vcf


class VcfFile(object):
    """Contains info about VCF file, and current position in file."""
    def __init__(
        self, file_name=None, file_=None, counts=None,
        reader=None, curr_rec=None, next_rec=None, curr_call=None,
        eof=False, has_curr_variant='absent', family_rel=None,
        curr_qc_values=None, is_rediscovery_file=False,
        sites_that_differ=None, unique_vars=None, absent_vars=None
    ):
        self.file_name = file_name
        self.file_ = file_
        self.counts = counts    # Category:count dict
        self.reader = reader    # PyVCF reader
        self.curr_rec = curr_rec    # Current record
        self.next_rec = next_rec    # Next record
        self.curr_call = curr_call   # Current call
        self.eof = eof    # End of file has been reached
        self.has_curr_variant = has_curr_variant  # Absent/low_qual/high_qual
        # DP, QUAL, etc. for current variant
        self.curr_qc_values = curr_qc_values
        # Is variant database used for rediscovery rate calculation
        self.is_rediscovery_file = is_rediscovery_file
        # Site: how site differs from others (string: string dict)
        self.sites_that_differ = sites_that_differ
        self.unique_vars = unique_vars  # Variants in this file not in others
        # Variants in other files not in this file
        self.absent_vars = absent_vars

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
    validate(args, parser)
    counts_per_file_categs = []  # Per-file statistics
    multi_sample_files = []  # Names of multi-sample VCF files
    vcf_file_objs = []
    empty_files = []

    with open(args.out_all_variants, 'w') as out_all:
        initialize_files(
            multi_sample_files, vcf_file_objs, empty_files, args,
            counts_per_file_categs, out_all
        )
        compare_variants(vcf_file_objs, args, out_all)

    for v in vcf_file_objs:
        v.file_.close()
        if not v.is_rediscovery_file:
            calculate_stats(v, vcf_file_objs)

    output_html_stats(
        counts_per_file_categs, vcf_file_objs, args, multi_sample_files,
        empty_files
    )


def parse_arguments(parser):
    """Add arguments to parser. Return parsed arguments."""
    parser.add_argument(
        'vcf_file_names', nargs='+',
        help='Names of VCF files to be processed (required)'
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
        '--out_all_variants', default='all_variants.txt',
        help='Name of output file containing all high-quality '
             'variants. (Default: all_variants.txt)'
    )
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()


def validate(args, parser):
    """Validate command line arguments."""
    message = ''
    for x in args.vcf_file_names + args.rediscovery_files:
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


def initialize_files(
    multi_sample_files, vcf_file_objs, empty_files, args,
    counts_per_file_categs, out_all
):
    """Create and initialize VcfFile objects and output file."""

    # Per-file statistics, in the order in which they will be displayed
    # in HTML table. (Categories starting with a capital letter are
    # displayed.) Stored in counts dict in VcfFile object.
    counts_per_file_categs += [
        'Total variants', 'Average DP', 'Average QUAL', 'average QD',
        'total_dp', 'count_dp', 'total_qual', 'count_qual', 'total_qd',
        'count_qd', 'High-quality variants', 'Homozygous reference', 'SNP',
        'Heterozygous', 'Homozygous alternate', 'Indel', 'Missing',
        'a>c', 'a>g', 'a>t', 'c>a', 'c>g', 'c>t', 'g>a', 'g>c', 'g>t', 't>a',
        't>c', 't>g', 'Transitions/ transversions',
        'Variants unique to this file',
        'Variants absent from this file present in other files'
    ]
    for f_name in args.rediscovery_files:
        counts_per_file_categs.append(
            'Rediscovery rate (non-indels): %s' % f_name
        )
        counts_per_file_categs.append('rediscovery_count_%s' % f_name)

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
                curr_qc_values={}, sites_that_differ={}, unique_vars=[],
                absent_vars=[]
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
    # First line of variant output file
    header = (
        'CHROM\tPOS\tREF\tALT\tGT\t%s\t%%Present\tAvg_QD\tAvg_DP\t'
        'Avg_QUAL\tAvg_AD0\tAvg_AD1\n' %
        ('\t'.join([v.file_name for v in vcf_file_objs]))
    )
    out_all.write(header)


def compare_variants(vcf_file_objs, args, out_all):
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
                parse_chromosome(v.next_rec.CHROM) ==
                parse_chromosome(min_chrom) and
                (not v.eof)
            )
        ])
        for v in vcf_file_objs:
            if (v.reader is None) or (v.next_rec is None) or (v.counts == {}):
                continue
            v.curr_rec = v.next_rec
            v.curr_call = None
            record = v.next_rec
            call = (
                record.genotype(v.reader.samples[0]) if (v.reader.samples)
                else None
            )
            if (
                parse_chromosome(record.CHROM) == parse_chromosome(min_chrom)
                and record.POS == min_pos and not v.eof
            ):
                if call:
                    assign_qc_values(v.curr_qc_values, record, call)
                    update_counts_all(v.counts, v.curr_qc_values)
                    is_high_qual = record.FILTER != ['LowQual']
                    v.curr_call = call

                    if is_high_qual:
                        v.has_curr_variant = 'high_qual'
                        update_counts_high_qual(v.counts, record, call)

                        for y in vcf_file_objs:
                            if (
                                y.is_rediscovery_file and
                                parse_chromosome(y.next_rec.CHROM) ==
                                parse_chromosome(min_chrom) and
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
        update_sites_that_differ(vcf_file_objs, min_chrom, min_pos)
        process_variant(vcf_file_objs, args, out_all)


def update_sites_that_differ(vcf_file_objs, chrom, pos):
    """Update dict of sites that differ between files."""
    if len(set(
        [v.has_curr_variant == 'absent' for v in vcf_file_objs]
    )) > 1:
        for y in vcf_file_objs:
            y.sites_that_differ['%s %d' % (chrom, pos)] = (
                'Variant %s<br />' % (
                    'present' if y.has_curr_variant != 'absent'
                    else 'absent'
                )
            )
    if len(set([
        (v.curr_rec.FILTER[0] if v.curr_rec.FILTER else 'NA')
        for v in vcf_file_objs if v.has_curr_variant != 'absent'
    ])) > 1:
        for y in vcf_file_objs:
            if y.has_curr_variant == 'absent':
                continue
            key = '%s %d' % (chrom, pos)
            value = (
                'FILTER: %s<br />(QUAL: %s)<br />' % (
                    'LowQual' if (y.curr_rec.FILTER == ['LowQual'])
                    else 'PASS',
                    ('%.2f' % y.curr_rec.QUAL) if y.curr_rec.QUAL
                    else 'None'
                )
            )
            if key in y.sites_that_differ:
                y.sites_that_differ[key] += '<br />%s' % value
            else:
                y.sites_that_differ[key] = value


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
            print 'excluding QUAL outlier: %d' % qc_values['qual']
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


def process_variant(vcf_file_objs, args, out_all):
    """Process a variant uniquely identified by CHROM and POS.
    Identify and process sub-variants."""
    files_per_variant = {}     # (REF, ALT, GT): list of VCF file objs

    for v in vcf_file_objs:
        if v.has_curr_variant != 'absent':
            call = v.curr_call
            record = v.curr_rec
            ref = str(record.REF)
            alt = str(record.ALT[0])
            gt = call['GT']

            if (ref, alt, gt) in files_per_variant:
                files_per_variant[(ref, alt, gt)].append(v)
            else:
                files_per_variant[(ref, alt, gt)] = [v]

    # Output each (REF, ALT, GT) combination at given position on chromosome
    for x in files_per_variant:
        output_ref_alt_gt(
            x, out_all, vcf_file_objs, args, files_per_variant
        )


def output_ref_alt_gt(
    ref_alt_gt, out_all, vcf_file_objs, args, files_per_variant
):
    """Process a variant uniquely identified by combination of
    CHROM, POS, REF, ALT, and GT. Write to file.
    """
    # List of VCF files that have given REF, ALT, GT
    ref_alt_gt_files = files_per_variant[ref_alt_gt]
    chrom = ref_alt_gt_files[0].curr_rec.CHROM
    pos = ref_alt_gt_files[0].curr_rec.POS
    ref, alt, gt = ref_alt_gt

    qc_averages = average_qc_values(ref_alt_gt_files)
    percent_present = divide(
        100 * len(ref_alt_gt_files),
        len([z for z in vcf_file_objs if not z.is_rediscovery_file]),
        -1.0
    )
    variant_in_file_str = ''
    for v in vcf_file_objs:
        if ref_alt_gt_files == [v]:
            v.counts['Variants unique to this file'] += 1
            v.unique_vars.append(
                '%s %s %s %s %s' % (chrom, str(pos), ref, alt, gt)
            )
        if v not in ref_alt_gt_files:
            v.counts[
                'Variants absent from this file present in other files'
            ] += 1
            v.absent_vars.append(
                '%s %s %s %s %s' % (chrom, str(pos), ref, alt, gt)
            )
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

    out_str = (
        '%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n'
        % (
            chrom, str(pos), ref, alt, gt, variant_in_file_str,
            percent_present, qc_averages['qd'], qc_averages['dp'],
            qc_averages['qual'], qc_averages['ad0'], qc_averages['ad1']
        )
    )
    out_all.write(out_str)


def divide(num, denom, default):
    """Return num/denom, or default if denom is 0."""
    if denom > 0:
        return float(num) / denom
    else:
        return default


def average_qc_values(vcf_file_objs):
    """Return dict containing average value for each QC category."""
    qc_avgs = {}
    qc_cats = [z for z in vcf_file_objs[0].curr_qc_values]
    for cat in qc_cats:
        qc_avgs[cat] = []
    for v in vcf_file_objs:
        # if v.has_curr_variant == 'high_qual':
        for cat in v.curr_qc_values:
            qc_avgs[cat].append(v.curr_qc_values[cat])

    for cat in qc_cats:
        if None in qc_avgs[cat]:
            qc_avgs[cat] = -1
        else:
            qc_avgs[cat] = divide(sum(qc_avgs[cat]), len(qc_avgs[cat]), -1)
    return qc_avgs


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
    out_f.write(table_title)
    out_f.write(
        '<table border="1" style="table-layout:fixed; width:'
        '%dpx; word-wrap:break-word"><tr>' %
        (175 * (len(dict_list) + 1))
    )
    for x in header_list:
        out_f.write('<th>%s</th>' % x)
    out_f.write('</tr>')
    for cat in cats:
        if 'chr' in cat:
            out_f.write('<tr><td style="text-align:center">%s</td>' % cat)
        else:
            out_f.write('<tr><th>%s</th>' % cat)
        for d in dict_list:
            out_f.write('<td>%s</td>' % d[cat])
        out_f.write('</tr>')
    out_f.write('</table>')


def format_counts(categories, counts):
    """Add normalized values to counts dict. Used for displaying
    HTML table.
    """
    if len(counts) == 0:
        return []
    counts_formatted = [dict(z) for z in counts]
    for cat in categories:
        for f in xrange(len(counts)):
            if (
                isinstance(counts[f][cat], int) or
                isinstance(counts[f][cat], float) and not
                isinstance(counts[f][cat], bool)
            ):
                if len(set([z[cat] for z in counts])) == 1:
                    if isinstance(counts[f][cat], float):
                        counts_formatted[f][cat] = '%.2f' % (counts[f][cat])
                    else:
                        counts_formatted[f][cat] = '%d' % (counts[f][cat])
                else:
                    if isinstance(counts[f][cat], float):
                        counts_formatted[f][cat] = (
                            '<font color = "blue">%.2f</font>' % (
                                counts[f][cat]
                            )
                        )
                    else:
                        counts_formatted[f][cat] = (
                            '<font color = "blue">%d</font>' % (
                                counts[f][cat]
                            )
                        )
    return counts_formatted


def ts_tv_ratio(counts):
    """Calculate ratio of transitions / transversions."""
    ts = ['a>g', 'g>a', 'c>t', 't>c']
    tv = ['a>c', 'c>a', 'a>t', 't>a', 'c>g', 'g>c', 'g>t',  't>g']
    return divide(
        sum([counts[x] for x in ts]), sum([counts[y] for y in tv]), -1.0
    )


def calculate_stats(vcf_file_obj, vcf_file_objs):
    """Calculate statistics for given file."""
    v = vcf_file_obj
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


def output_html_stats(
    counts_per_file_categs, vcf_file_objs, args, multi_sample_files,
    empty_files
):
    """Output per-file statistics to HTML file."""
    counts_formatted = format_counts(
        counts_per_file_categs,
        [v.counts for v in vcf_file_objs if not v.is_rediscovery_file]
    )
    categs_pre_filter = ['Total variants', 'Average DP', 'Average QUAL']

    # Per unique combination of CHROM and POS
    categs_post_filter_1 = [
        'High-quality variants', 'Homozygous reference', 'SNP',
        'Heterozygous', 'Homozygous alternate', 'Indel', 'Missing',
        'Transitions/ transversions'
    ]
    for f_name in args.rediscovery_files:
        categs_post_filter_1.append(
            'Rediscovery rate (non-indels): %s' % f_name
        )
    # Per unique combination of CHROM, POS, REF, ALT, GT
    categs_post_filter_2 = [
        'Variants unique to this file',
        'Variants absent from this file present in other files',
    ]

    page_title = 'VCF Comparison'
    with open(args.out_stats, 'w') as out_f:
        out_f.write(
            '<!DOCTYPE html><html><head><title>%s</title></head><body>'
            '<h1>%s</h1><br />' % (page_title, page_title)
        )
        print_dict_list(
            out_f, '<h2>Statistics per file, all variants,<br />per unique '
            'combination of CHROM and POS:</h2>', counts_formatted,
            categs_pre_filter, [''] + [
                v.file_name for v in vcf_file_objs if not
                v.is_rediscovery_file
            ]
        )
        if any([v.sites_that_differ for v in vcf_file_objs]):
            print_dict_list(
                out_f, '<h3>Sites that differ between files:</h3>',
                [z.sites_that_differ for z in vcf_file_objs],
                order_sites([z for z in v.sites_that_differ]),
                ['Site'] + [
                    v.file_name for v in vcf_file_objs if not
                    v.is_rediscovery_file
                ]
            )
        print_dict_list(
            out_f, '<br /><h2>Statistics per file, high-quality variants '
            'only,<br />per unique combination of CHROM and POS:</h2>',
            counts_formatted, categs_post_filter_1,
            [''] + [
                v.file_name for v in vcf_file_objs if not
                v.is_rediscovery_file
            ]
        )
        print_dict_list(
            out_f, '<br /><h2>Statistics per file, high-quality variants '
            'only,<br />per unique combination of CHROM, POS, REF, ALT, '
            'and GT:</h2>', counts_formatted, categs_post_filter_2,
            [''] + [
                v.file_name for v in vcf_file_objs if not
                v.is_rediscovery_file
            ]
        )
        if len(multi_sample_files) > 0:
            out_f.write(
                '<p><b>Multi-sample files excluded from ' +
                'analysis: </b>%s</p>' % (', '.join(multi_sample_files)))
        if len(empty_files) > 0:
            out_f.write(
                '<p><b>Empty/ invalid files excluded from ' +
                'analysis: </b>%s</p>' % (', '.join(empty_files)))
        for v in vcf_file_objs:
            if v.unique_vars:
                out_f.write(
                    '<p><b>Variants unique to %s (CHROM, POS, REF, ALT, GT):'
                    '<br /></b> %s</p>' % (
                        v.file_name, '<br />'.join(v.unique_vars)
                    )
                )
            if v.absent_vars:
                out_f.write(
                    '<p><b>Variants absent from %s (CHROM, POS, REF, ALT, '
                    'GT):<br /></b> %s</p>' % (
                        v.file_name, '<br />'.join(v.absent_vars)
                    )
                )
        out_f.write('</body></html>')


def order_sites(list_of_sites):
    """Order list of sites [chrV n, chrY p, ...], where chrV and chrY are
    CHROM values, and n and p are POS values.
    """
    l_ = list(list_of_sites)
    l_ordered = []
    while len(l_) > 0:
        min_chrom = min_chromosome(z.split(' ')[0] for z in l_)
        for y in sorted([
            int(z.split(' ')[1]) for z in l_ if z.split(' ')[0] == min_chrom
        ]):
            site = '%s %d' % (min_chrom, y)
            l_ordered.append(site)
            l_.remove(site)
    return l_ordered


if __name__ == '__main__':
    sys.exit(main())
