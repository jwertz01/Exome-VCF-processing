"""Compares VCF files obtained from a parent-offspring trio.
Produces statistics and variant comparison files, with analysis
of inheritance."""

import os
import sys
import argparse
import matplotlib.pyplot as plt
import vcf


class VcfFile(object):
    def __init__(
        self=None, file_name=None, file_=None, counts=None,
        reader=None, curr_rec=None, next_rec=None, curr_call=None,
        eof=False, has_curr_variant='absent', family_rel=None
    ):
        self.file_name = file_name
        self.file_ = file_
        self.counts = counts    # category:count map
        self.reader = reader    # PyVCF reader
        self.curr_rec = curr_rec    # current record
        self.next_rec = next_rec    # next record
        self.curr_call = curr_call   # current call
        self.eof = eof    # end of file has been reached
        self.has_curr_variant = has_curr_variant   # absent/low_qual/high_qual
        self.family_rel = family_rel  # parent1/parent2/child

    def __str__(self):
        return (
            'Filename: %s\nCounts: %s\nCurr rec: %s\nNext rec: %s\nEOF: %s'
            '\nHas curr variant: %s\nFamily rel: %s\n' %
            (
                self.file_name, self.counts, self.curr_rec, self.next_rec,
                self.eof, self.has_curr_variant, self.family_rel
            )
        )


def main():
    args = add_arguments(argparse.ArgumentParser())
    if (
        args.parent_1_file_names is None or
        args.parent_2_file_names is None or
        args.child_file_names is None
    ):
        raise ValueError(
            'Parent 1, Parent 2, and child file names must be specified.'
        )
    multi_sample_files = []  # names of multi-sample VCF files
    vcf_file_objs = []

    # qc plot x-axis intervals
    qc_intervals = {
        'dp': xrange(0, 55, 5), 'qd': xrange(0, 20, 2),
        'qual': xrange(0, 260, 26)
    }

    # categories displayed in qc plots
    qc_categories = [
        'All', 'De novo', 'Inherited', 'Invalid', 'Internal conflict'
    ]
    qc_counts = {}
    # categories displayed in flowchart in html file
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

    # per-file statistics, in the order in which they will be displayed
    # in html table
    categories = [
        'Total variants', 'High-quality variants',
        'Homozygous reference', 'SNP', 'Heterozygous', 'Homozygous alternate',
        'Indel', 'Missing', 'A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T', 'G>A',
        'G>C', 'G>T', 'T>A', 'T>C', 'T>G'
    ]

    # initialize readers, counts
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
                file_=in_f, counts={}, reader=reader, next_rec=reader.next()
            )
            if f_name in args.parent_1_file_names:
                v.family_rel = 'Parent 1'
            elif f_name in args.parent_2_file_names:
                v.family_rel = 'Parent 2'
            else:
                v.family_rel = 'Child'
            for cat in categories:
                v.counts[cat] = 0
            v.file_name = '%s (%s)' % (f_base_name, v.family_rel)
            vcf_file_objs.append(v)

    if args.create_qc_plots:
        for i in qc_intervals:
            qc_counts[i] = {}
            for j in qc_categories:
                qc_counts[i][j] = {}
                for k in qc_intervals[i]:
                    qc_counts[i][j][k] = 0
    for x in overall_counts_categs:
        overall_counts[x] = 0

    # iterate over vcf files, write variant comparison output files,
    # fill in counts
    header = (
        'CHROM\tPOS\tREF\tALT\tGT\tInheritance_category\t'
        'Inheritance_subcategory\t%s\t%%Present\tAve_QD'
        '\tAve_DP\tAve_QUAL\tAve_AD0\tAve_AD1\n' %
        ('\t'.join([v.file_name for v in vcf_file_objs]))
    )
    with open(args.out_all_variants, 'w') as out_all:
        out_all.write(header)
        compare_variants(
            vcf_file_objs, qc_counts, overall_counts, args, out_all
        )

    for v in vcf_file_objs:
        v.file_.close()

    # output qc plots
    if args.create_qc_plots:
        output_qc_plot(
            qc_categories, 'Read depth', qc_counts['dp'], args.depth_plot,
            qc_intervals['dp']
        )
        output_qc_plot(
            qc_categories, 'Quality by depth', qc_counts['qd'], args.qd_plot,
            qc_intervals['qd']
        )
        output_qc_plot(
            qc_categories, 'Quality', qc_counts['qual'], args.qual_plot,
            qc_intervals['qual']
        )

    # output html stats file
    counts_norm = normalize_counts(
        categories, [v.counts for v in vcf_file_objs]
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
            out_f, 'Counts Per File', counts_norm, categories,
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


def add_arguments(parser):
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
        '--out_all_variants', default='out_all_variants.txt',
        help='Name of output file containing all high-quality '
             'variants. (Default: all_pass_qc.txt)'
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
    return parser.parse_args()


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
            qd = record.INFO['QD']
            dp = call['DP']
            qual = record.QUAL
            if (
                record.CHROM == min_chrom and record.POS == min_pos
                and not v.eof
            ):
                v.counts['Total variants'] += 1
                is_high_qual = (
                    qd >= args.min_qd and dp >= args.min_dp and
                    qual >= args.min_qual
                )
                if is_high_qual:
                    v.has_curr_variant = 'high_qual'
                    update_counts(v.counts, v.reader, v.next_rec)
                else:
                    v.has_curr_variant = 'low_qual'
                v.curr_call = call
                try:
                    v.next_rec = v.reader.next()
                except StopIteration:
                    v.eof = True
            else:
                v.has_curr_variant = 'absent'

        # write single variant (at min_chrom and min_pos) to output files,
        # update counts
        process_variant(
            vcf_file_objs, qc_counts, overall_counts, args, out_all
        )


def process_variant(
    vcf_file_objs, qc_counts, overall_counts, args, out_all
):
    """Output a single variant to output files, and update counts."""
    if all(
        v.has_curr_variant in ['low_qual', 'absent'] for v in
        vcf_file_objs
    ):
        return

    files_per_variant = {}     # (ref, alt, gt): list of VCF file objs
    for v in vcf_file_objs:
        if v.has_curr_variant != 'absent':
            ref = str(v.curr_rec.REF)
            alt = str(v.curr_rec.ALT[0])
            gt = v.curr_call['GT']
            if (ref, alt, gt) in files_per_variant:
                files_per_variant[(ref, alt, gt)].append(v)
            else:
                files_per_variant[(ref, alt, gt)] = [v]

    # output each (ref, alt, gt) combination at given position on chromosome
    for x in files_per_variant:
        output_ref_alt_gt(
            vcf_file_objs, qc_counts, overall_counts, args, out_all, x,
            files_per_variant
        )


def output_ref_alt_gt(
    vcf_file_objs, qc_counts, overall_counts, args, out_all, ref_alt_gt,
    files_per_variant
):
    """Process a variant uniquely identified by combination of
    CHROM, POS, REF, ALT, and GT): Write to file, and update
    overall_counts.
    """
    # list of VCF files that have given ref, alt, gt
    if all([v.has_curr_variant != 'high_qual' for v in vcf_file_objs]):
        return
    ref_alt_gt_files = files_per_variant[ref_alt_gt]
    chrom = ref_alt_gt_files[0].curr_rec.CHROM
    pos = ref_alt_gt_files[0].curr_rec.POS
    ref, alt, gt = ref_alt_gt

    qc_averages = {'qd': 0, 'dp': 0, 'qual': 0, 'ad0': 0, 'ad1': 0}
    count_high_qual = 0
    for v in ref_alt_gt_files:
        if v.has_curr_variant == 'high_qual':
            call = v.curr_call
            record = v.curr_rec
            qc_averages['qd'] += record.INFO['QD']
            qc_averages['dp'] += call['DP']
            qc_averages['qual'] += record.QUAL
            qc_averages['ad0'] += call['AD'][0]
            qc_averages['ad1'] += call['AD'][1]
            count_high_qual += 1
    for x in qc_averages:
        if count_high_qual == 0:
            qc_averages[x] = -1
        else:
            qc_averages[x] /= float(count_high_qual)

    percent_present = 100.0 * len(ref_alt_gt_files) / len(vcf_file_objs)

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

    overall_counts['Pass QC'] += 1

    # output indels
    if True in [
        v.curr_rec.is_indel for v in vcf_file_objs if
        v.has_curr_variant != 'absent'
    ]:
        out_all.write(
            '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t'
            '%.2f\t%.2f\n' %
            (
                chrom, str(pos), ref, alt, gt, '[Indel]', 'NA',
                variant_in_file_str, percent_present, qc_averages['qd'],
                qc_averages['dp'], qc_averages['qual'], qc_averages['ad0'],
                qc_averages['ad1']
            )
        )
        overall_counts['Indel'] += 1
        return
    overall_counts['No indel'] += 1

    # find internal disagreements within parent 1 files/
    # parent 2 files/ child files
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
    if conflict_exists:
        cat = 'Internal conflict'
        subcat = 'NA'
    else:
        cat, subcat = inheritance(
            *[example_files[family_rel] for family_rel in family_rels], pos=pos
        )
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
    overall_counts[cat] += 1
    if cat == 'Inherited':
        overall_counts[subcat] += 1

    # update plot data
    if args.create_qc_plots:
        update_qc_counts(
            ave_qc_vals(ref_alt_gt_files), qc_counts, len(files_per_variant),
            cat
        )


def internal_conflict_exists(vcf_file_objs, family_rel, files_per_variant):
    """Returns whether there is an internal conflict within files with
    given family_rel."""
    conflict_exists = False
    first_instance = None
    for x in files_per_variant:  # x == ref, alt, gt
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


def update_qc_counts(ave_qc_vals, qc_counts, count_ref_alt_objs, categ):
    """Update qc_counts dict."""
    for qc_cat in ave_qc_vals:
        intervals = qc_counts[qc_cat]['All']
        # normalize by number of distinct ref/alt/gt combinations
        qc_val = ave_qc_vals[qc_cat] / count_ref_alt_objs
        for x in intervals:
            if qc_val >= x:
                qc_counts[qc_cat][categ][x] += 1


def ave_qc_vals(vcf_file_objs):
    """Returns dict containing average value for each QC category,
    among files with given family_rel (child/ parent1/ parent2)."""
    qc_vals = {}
    qc_cats = ['dp', 'qd', 'qual']
    for cat in qc_cats:
        qc_vals[cat] = []
    for v in vcf_file_objs:
        qc_vals['dp'].append(v.curr_call['DP'])
        qc_vals['qd'].append(v.curr_rec.INFO['QD'])
        qc_vals['qual'].append(v.curr_rec.QUAL)
    for cat in qc_cats:
        qc_vals[cat] = (
            (float(sum(qc_vals[cat])) / len(qc_vals[cat])) if qc_vals[cat]
            else -1
        )
    return qc_vals


def min_chromosome(chroms):
    """Args: Chromosome names ['chrN', 'chrP', ...]. Returns name of
    chromosome which would appear first in ordered file.
    """
    chroms_processed = list(chroms)
    for i in xrange(len(chroms_processed)):
        chroms_processed[i] = (
            chroms_processed[i][chroms_processed[i].index('chr')+3:])
        try:
            chroms_processed[i] = int(chroms_processed[i])
        except ValueError:
            pass
    return chroms[chroms_processed.index(min(chroms_processed))]


def update_counts(counts, reader, record):
    """Update counts map given a VCF record. Assumes high-quality
    variant.

    Args:
        counts (dict): Dict of category:count.
        reader (VCF reader): VCF reader for file.
        record (VCF reader record): Specific VCF record.
    """
    call = record.genotype(reader.samples[0])
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


def print_dict_list(
    out_f, table_title, dict_list, cats, header_list, table_message
):
    """Create HTML file and print a list of category:value
    dictionaries as an HTML table.

    Args:
        out_f (str): File to write to.
        table_title (str): HTML table title.
        dict_list (list of dict): Category:value dicts.
        cats (list of string): Ordered categories for stats output.
        header_list (str): List to group categories by.
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
    # normalized table entries at least flag_dist
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
                    (not isinstance(counts[f][cat], int)) or
                    (isinstance(counts[f][cat], bool)) or
                    (min_val == 0)
                ):
                    continue
                counts_norm[f][cat] = (
                    counts[f][cat] / float(min_val) *
                    (
                        float(counts[min_index]['High-quality variants']) /
                        counts[f]['High-quality variants']
                    )
                )
                if abs(1.0 - counts_norm[f][cat]) < flag_dist:
                    counts_norm[f][cat] = (
                        '%d (%.2f)' % (counts[f][cat], counts_norm[f][cat])
                    )
                else:
                    counts_norm[f][cat] = (
                        '%d (<font color = "blue">%.2f</font>)'
                        % (counts[f][cat], counts_norm[f][cat])
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


def inheritance(parent1_file, parent2_file, child_file, pos):
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

    # process variants for which there is a ref or alt conflict between
    # parent1 and parent2, or parent1 and child, or parent2 and child
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


if __name__ == '__main__':
    main()
