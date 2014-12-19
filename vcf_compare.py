"""Compares VCF files. Produces statistics and variant comparison files."""

import os
import sys
import argparse
import vcf
import matplotlib.pyplot as plt

class VcfFile:
	def __init__(self, file_name, file, reader, next_rec):
		self.file_name = file_name
		self.file = file
		self.counts = {}	# category:count map
		self.reader = reader	# PyVCF reader
		self.curr_rec = None	# previous record
		self.next_rec = next_rec	# current record
		self.curr_call = None	# current call
		self.eof = False	# end of file has been reached
		self.has_curr_variant = 'absent'	# absent/low_qual/high_qual
	
	def __str__(self):
		return (
			'Filename: %s\nCounts: %s\nCurr rec: %s\nNext rec: %s\nEOF: %s' 
			'\nHas curr variant: %s\nFamily rel: %s\n' % 
			(self.file_name, self.counts, self.curr_rec, self.next_rec, 
			self.eof, self.has_curr_variant, self.family_rel))
		
def main():
	parser = argparse.ArgumentParser()	
	parser.add_argument('vcf_file_names', nargs = '+')
	parser.add_argument('--create_qc_plots', action = 'store_true')
	parser.add_argument('--out_stats', default = 'vcf_stats.html')
	parser.add_argument('--depth_plot', default = 'dp_plot.png')
	parser.add_argument('--qd_plot', default = 'qd_plot.png')
	parser.add_argument('--qual_plot', default = 'qual_plot.png')
	
	# containing all high-quality variants
	parser.add_argument('--out_all_variants', default = 'all_variants.txt')
			
	# containing variants such that at least one file has high-quality 
	# variant and another file does not have variant at chrom and pos.
	parser.add_argument('--out_absent', default = 'conflicts_absent.txt')
	
	# containing variants such that at least one file has high-quality
	# variant and another file has different high-quality ref/alt combination 
	# at chrom and pos.
	parser.add_argument(
		'--out_multiple_vars', default = 'conflicts_multiple_vars.txt')
	parser.add_argument(
		'--out_multiple_vars_no_indel', 
		default = 'conflicts_multiple_vars_no_indel.txt')
	
	# min quality by depth, quality, read depth for inclusion 
	# in above variant output files
	parser.add_argument('--min_qd', type = int, default = 3) 
	parser.add_argument('--min_qual', type = int, default = 80) 
	parser.add_argument('--min_dp', type = int, default = 15)

	args = parser.parse_args()
	multi_sample_files = []	#names of multi-sample VCF files
	vcf_file_objs = []
	
	#qc plot x-axis intervals
	qc_intervals = {
		'dp': xrange(0, 55, 5), 'qd': xrange(0, 20, 2), 
		'qual': xrange(0, 260, 26)}
	qc_categories = [
		'all', 'het_all_agree', 'het_disagree', 'hom_all_agree', 
		'hom_disagree']
	qc_counts = {}
	overall_counts_categs = ['Pass_QC', 'All_agree', 'Conflicts', 
	'Conflicts_absent', 'Conflicts_multiple', 'Conflicts_multiple_indel',
	'Conflicts_multiple_no_indel', 'Low_QC', 'Agree_GT_disagree_ALT',
	'Disagree_GT_disagree_ALT', 'Disagree_GT_agree_ALT', 'Predict_HET', 
	'Predict_HOM_ALT', 'No_prediction']
	overall_counts = {}
	
	#per-file statistics, in the order in which they will be displayed
	#in html table
	categories=[
		'Total variants', 'High-quality variants',
		'Homozygous reference', 'SNP', 'Heterozygous', 'Homozygous alternate',
		'Indel', 'Missing', 'A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T', 'G>A', 
		'G>C', 'G>T', 'T>A', 'T>C', 'T>G',
		'Variants absent from this file present in other files',
		'High-quality variants present that conflict with ' 
		'variants present in other files']
			
	#initialize readers, counts
	for f_name in args.vcf_file_names:
		in_f = open(f_name)
		reader = vcf.Reader(in_f)
		f_base_name = os.path.basename(f_name)		
		if len(reader.samples) > 1:
			multi_sample_files.append(f_base_name)
			in_f.close()
		else:
			v = VcfFile(f_base_name, in_f, reader, reader.next())
			for cat in categories:
				v.counts[cat] = 0
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
			
	#iterate over vcf files, write variant comparison output files, 
	#fill in counts	
	header = (
		'CHROM\tPOS\tREF\tALT\tGT\t%s\t%%Present\tAve_QD\tAve_DP\tAve_AD0'
		'(high-qual)\tAve_AD1(high-qual)\tAve_QUAL' %
		('\t'.join([v.file_name for v in vcf_file_objs])))
	with open(args.out_all_variants, 'w') as out_all, \
	open(args.out_absent, 'w') as out_absent, \
	open(args.out_multiple_vars, 'w') as out_multiple, \
	open(args.out_multiple_vars_no_indel, 'w') as out_multiple_no_indel:
		out_all.write(header + '\n')		
		out_absent.write(header + '\n')
		out_multiple.write(header + '\n')
		out_multiple_no_indel.write(header + '\tSummary\tConclusion\n')
	
		compare_variants(
			vcf_file_objs, qc_counts, overall_counts, args, 
			out_all, out_absent, out_multiple, out_multiple_no_indel)
		
	for v in vcf_file_objs:
		v.file.close()
	
	#output qc plots
	if args.create_qc_plots:
		output_qc_plot(
			qc_categories, 'Read depth', qc_counts['dp'], args.depth_plot, 
			qc_intervals['dp'])
		output_qc_plot(
			qc_categories, 'Quality by depth', qc_counts['qd'], args.qd_plot, 
			qc_intervals['qd'])
		output_qc_plot(
			qc_categories, 'Quality', qc_counts['qual'], args.qual_plot, 
			qc_intervals['qual'])
		
 	#output html stats file
 	counts_norm = normalize_counts(
 		categories, [v.counts for v in vcf_file_objs])
	with open(args.out_stats, 'w') as out_f:
		page_title = 'VCF Comparison'
		out_f.write(
			'<!DOCTYPE html><html><head><title>%s</title></head><body>'
			'<h1>%s</h1>' % (page_title, page_title))
		table_message = (
			'<p>The bottom two rows are per unique combination of CHROM, '
			'POS, REF, ALT, and GT. The rest are per unique combination '
			'of CHROM and POS. </p>'
			'<p>Normalized with respect to the number of '
			'high-quality variants, and to the value in the leftmost '
			'numeric column. Normalized values less than 0.95 or '
			'greater than 1.05 are shown in blue. Counts shown (excluding '
			'total variants) are after quality filtering.</p>')
		print_dict_list(
			out_f, 'Counts Per File', counts_norm, categories, 
			[v.file_name for v in vcf_file_objs], table_message)
		if len(multi_sample_files) > 0:
			out_f.write(
				'<p><b>Multi-sample files excluded from ' +
				'analysis: </b>%s</p>' % (', '.join(multi_sample_files)))
		out_f.write(
			'<h2>Overall Counts Per Unique Combination of CHROM, POS, REF, '
			'ALT, and GT</h2>')
		out_f.write(
			'<p>Each category (except Pass_QC) is a subcategory of the '
			'rightmost category directly above it.</p>')
		out_f.write('<table><tr>')
		for x in overall_counts_categs:
			if (
				x in ['All_agree', 'Conflicts_absent', 
				'Conflicts_multiple_indel', 'Low_QC', 'Predict_HET']):
				out_f.write('</tr><tr>')
			out_f.write(
				'<td style=padding-right:2em><b>%s</b>: %d</td>' % 
				(x, overall_counts[x]))
		out_f.write('</tr></table>')
		out_f.write('</body></html>')
		

def compare_variants(
	vcf_file_objs, qc_counts, overall_counts, args, out_all, 
	out_absent, out_multiple, out_multiple_no_indel):
	"""Iterate through all VCF files simultaneously. Write per-file 
	variant info to output files and update overall statistics (counts). 
	"""
	for v in vcf_file_objs:
		if (v.reader is None) or (v.next_rec is None) or (v.counts == {}):
			return
			
	min_chrom = min_chromosome(
		[v.next_rec.CHROM for v in vcf_file_objs if not v.eof])
	print 'Processing %s...' % min_chrom
	while not (all([v.eof for v in vcf_file_objs])):
		new_min_chrom = min_chromosome(
			[v.next_rec.CHROM for v in vcf_file_objs if not v.eof])
		if new_min_chrom != min_chrom:
			print 'Processing %s...' % new_min_chrom
			min_chrom = new_min_chrom
		min_pos = min(
			[v.next_rec.POS for v in vcf_file_objs if 
			v.next_rec.CHROM == min_chrom and not v.eof])
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
				and not v.eof):
				v.counts['Total variants'] += 1		
				is_high_qual = (
					qd >= args.min_qd and dp >= args.min_dp and 
					qual >= args.min_qual)
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
				
		#write single variant (at min_chrom and min_pos) to output files,
		#update counts
		process_variant(
			vcf_file_objs, out_all, out_absent, out_multiple, 
			out_multiple_no_indel, qc_counts, overall_counts, args)
	
#used with process_variant function		
class RefAlt:
	"""Object representing a variant (uniquely identified by combination of
	CHROM, POS, REF, ALT, and GT).
	"""
	def __init__(self, num_records):
		self.total_qd = 0
		self.total_qual = 0
		self.total_dp = 0
		self.total_ad0 = 0
		self.total_ad1 = 0
		self.in_file = ['absent'] * num_records	
		
def process_variant(
	vcf_file_objs, out_all, out_absent, out_multiple, out_multiple_no_indel, 
	qc_counts, overall_counts, args):
	"""Output a single variant to output files, and update counts.
	"""
	high_qual_present = (
		not all(v.has_curr_variant in ['low_qual', 'absent'] for v in 
		vcf_file_objs))
	
	#find each (ref, alt) combination at given position
	#on chromosome
	dict_var_stats = {} 	#(ref, alt) : RefAlt obj

	for i, v in enumerate(vcf_file_objs):
		if v.has_curr_variant == 'absent':
			continue		
		call = v.curr_call
		record = v.curr_rec
		ref = str(record.REF)
		alt = str(record.ALT[0])
		dp = call['DP']
		qd = record.INFO['QD']
		qual = record.QUAL
		gt = call['GT']
		ad0 = call['AD'][0]
		ad1 = call['AD'][1]
		#update dict_var_stats
		if high_qual_present:			
			if not (ref, alt, gt) in dict_var_stats:
				dict_var_stats[(ref, alt, gt)] = RefAlt(len(vcf_file_objs))
			ref_alt_obj = dict_var_stats[(ref, alt, gt)]
			dict_var_stats[(ref, alt, gt)].in_file[i] = v.has_curr_variant					
			if v.has_curr_variant == 'high_qual':
				ref_alt_obj.total_qd += qd
				ref_alt_obj.total_dp += dp
				ref_alt_obj.total_ad0 += ad0
				ref_alt_obj.total_ad1 += ad1
				ref_alt_obj.total_qual += qual
	multiple_vars = (
		len([z for z in dict_var_stats if 'high_qual' in 
		dict_var_stats[z].in_file]) > 1)
		
	#update qc counts
	if args.create_qc_plots:
		all_agree = (
			all(v.has_curr_variant != 'absent' for v in vcf_file_objs) and 
			(not multiple_vars))
		for v in vcf_file_objs:
			if v.has_curr_variant == 'absent':
				continue
			call = v.curr_call
			record = v.curr_rec
			try:
				dp = call['DP']
			except TypeError:
				dp = None 
			qd = record.INFO['QD']
			qual = record.QUAL
			update_qc_counts(record, all_agree, dp, qc_counts['dp'])
			update_qc_counts(record, all_agree, qd, qc_counts['qd'])
			update_qc_counts(record, all_agree, qual, qc_counts['qual'])
						
	if high_qual_present:		
		#output each (ref, alt) combination at given position on chromosome
		for x in dict_var_stats:
			 output_ref_alt(
			 	x, out_all, out_absent, out_multiple, out_multiple_no_indel,
			 	dict_var_stats, vcf_file_objs, overall_counts)
	
	
def update_qc_counts(record, all_agree, qc_val, qc_counts):
	"""Update qc_counts dict.
	"""
	intervals = [z for z in qc_counts['all']]
	for x in intervals:
		if (qc_val is not None) and qc_val >= x:
			qc_counts['all'][x] += 1
			if all_agree:
				qc_counts['hom_all_agree'][x] += record.num_hom_alt
				qc_counts['het_all_agree'][x] += record.num_het
			else:
				qc_counts['hom_disagree'][x] += record.num_hom_alt
				qc_counts['het_disagree'][x] += record.num_het	
			
			
def output_ref_alt(
	ref_alt, out_all, out_absent, out_multiple, out_multiple_no_indel,
	dict_var_stats, vcf_file_objs, overall_counts):
	"""Process a variant uniquely identified by combination of
	CHROM, POS, REF, ALT, and GT): Write to file, and update
	overall_counts.
	"""
	ref_alt_obj = dict_var_stats[ref_alt]
	if all (z != 'high_qual' for z in ref_alt_obj.in_file):
		return
		
	i = [v.has_curr_variant for v in vcf_file_objs].index('high_qual')
	chrom = vcf_file_objs[i].curr_rec.CHROM
	pos = vcf_file_objs[i].curr_rec.POS
	
	ref = ref_alt[0]
	alt = ref_alt[1]
	gt = ref_alt[2]
	count_high_qual = ref_alt_obj.in_file.count('high_qual')
	ave_qd = float(ref_alt_obj.total_qd) / count_high_qual
	ave_dp = float(ref_alt_obj.total_dp) / count_high_qual
	ave_ad0 = float(ref_alt_obj.total_ad0) / count_high_qual
	ave_ad1 = float(ref_alt_obj.total_ad1) / count_high_qual
	ave_qual = float(ref_alt_obj.total_qual) / count_high_qual
	multiple_variants = len([z for z in dict_var_stats if 'high_qual' in 
		dict_var_stats[z].in_file]) > 1
	percent_present = (
		100.0 * (ref_alt_obj.in_file.count('low_qual') + 
		ref_alt_obj.in_file.count('high_qual')) / len(ref_alt_obj.in_file))
	variant_in_file_str = '\t'.join(ref_alt_obj.in_file)
	for i, in_file in enumerate(ref_alt_obj.in_file):
		if in_file == 'absent':
			#same chrom, pos but different ref/ alt counts as absent.
			vcf_file_objs[i].counts['Variants absent from this file '
			'present in other files'] += 1
		if (
			in_file == 'high_qual' and percent_present != 100
			and multiple_variants):
			vcf_file_objs[i].counts['High-quality variants present that '
			'conflict with variants present in other files'] += 1
	out_str = (
		'%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' 
		% (chrom, str(pos), ref, alt, gt, variant_in_file_str, 
		percent_present, ave_qd, ave_dp, ave_ad0, ave_ad1, ave_qual))
	overall_counts['Pass_QC'] += 1
	out_all.write(out_str)	
	
	if percent_present == 100:
		overall_counts['All_agree'] += 1
	else:
		overall_counts['Conflicts'] += 1
		if multiple_variants: 
			overall_counts['Conflicts_multiple'] += 1
			out_multiple.write(out_str)
			if all([len(r) == len(a) for (r, a, g) in dict_var_stats]):
				overall_counts['Conflicts_multiple_no_indel'] += 1
				out_multiple_no_indel.write(out_str[:-1])
				write_multiple_vars(
					dict_var_stats, out_multiple_no_indel, overall_counts)
			else:
				overall_counts['Conflicts_multiple_indel'] += 1
		else:	
			overall_counts['Conflicts_absent'] += 1
			out_absent.write(out_str)


def write_multiple_vars(dict_var_stats, out_multiple, overall_counts):
	"""Write variant info (summary and conclusion) to out_multiple file,
	and update overall_counts dict.
	"""
	#write conflict categories
	conflicts = []
	if len(set([r for (r, a, g) in dict_var_stats])) > 1:
		conflicts.append('REF')
	if len(set([a for (r, a, g) in dict_var_stats])) > 1:
		conflicts.append('ALT')
	if len(set([g for (r, a, g) in dict_var_stats])) > 1:
		conflicts.append('GT')

	#write most frequent values
	plur = []
	if 'REF' in conflicts:
		common_ref = plurality(dict_var_stats, 0)
		plur += common_ref
	if 'ALT' in conflicts:
		common_alt = plurality(dict_var_stats, 1)
		plur += common_alt
	if 'GT' in conflicts:
		common_gt = plurality(dict_var_stats, 2)
		plur += common_gt
	out_multiple.write(
		'\t%s conflict.' % ', '.join(conflicts) + ((' Most frequent: %s.' % 
		', '.join(plur)) if (plur != []) else ''))
		
	#write AD ratio summary
	if 'GT' in conflicts:
		ad_ratio = (
			float(sum([dict_var_stats[x].total_ad0 for x in dict_var_stats])) 
			/ sum([dict_var_stats[y].total_ad1 for y in dict_var_stats]))
		ad_ratio_categ = 'other'
		if ad_ratio <= .1:
			ad_ratio_categ = 'HOM-ALT'
		elif ad_ratio > 0.8 and ad_ratio < (1/.8):
			ad_ratio_categ = 'HET'
		elif (ad_ratio > 10):
			ad_ratio_categ = 'HOM-REF'
		if ad_ratio_categ != 'other':	
			out_multiple.write(
				' AD0/AD1 ratio indicates %s.' % ad_ratio_categ)
		
	#write QC summary
	max_qd_ref_alts = max(
		(((float(dict_var_stats[z].total_qd) / dict_var_stats[
		z].in_file.count('high_qual')) if (dict_var_stats[z
		].in_file.count('high_qual') > 0) else -1) for z in dict_var_stats))

	if max_qd_ref_alts < 5:
		out_multiple.write(' Low QD.')
	elif max_qd_ref_alts < 10:
		out_multiple.write(' Borderline QD.')

	#write conclusion
	if (max_qd_ref_alts < 10):
		overall_counts['Low_QC'] += 1
		out_multiple.write('\tNeeds validation (low/borderline QC).')
	elif 'REF' in conflicts:
		out_multiple.write('\tREF conflict(?)')
	elif 'GT' in conflicts and ad_ratio_categ == 'HOM-REF':
		out_multiple.write('\tHOM-REF AD ratio(?)')
	elif 'GT' in conflicts and 'ALT' in conflicts:
		overall_counts['Disagree_GT_disagree_ALT'] += 1
		out_multiple.write('\tDisagree GT, ALT.')
	elif 'GT' in conflicts:
		overall_counts['Disagree_GT_agree_ALT'] += 1
		out_multiple.write('\tDisagree GT, agree ALT.')
		if len(common_gt) == 1 and ((common_gt[0] == '0/1' and ad_ratio_categ 
		== 'HET') or (common_gt[0] == '1/1' and	ad_ratio_categ == 'HOM-ALT')):
			if ad_ratio_categ == 'HOM-ALT':
				overall_counts['Predict_HOM_ALT'] += 1
			elif ad_ratio_categ == 'HET':
				overall_counts['Predict_HET'] += 1
			out_multiple.write(' Predict %s.' % ad_ratio_categ)
		else:
			overall_counts['No_prediction'] += 1
	elif 'ALT' in conflicts:
		overall_counts['Agree_GT_disagree_ALT'] += 1
		out_multiple.write('\tAgree GT, disagree ALT.')
	else:
		out_multiple.write('\tNo conflicts(?)')
	out_multiple.write('\n')
	

def plurality(tuple_dict, n):
	""" 
		Args:
			tuple_dict (dict of (a, b, c) : RefAlt obj): Has in_file list
			n (int): Position in tuple_dict key
		Returns:
			Element(s) of (a, b, c)[n] that occurs most frequently in files
	"""
	l = (
		[[(a, b, c)[n]] * (tuple_dict[(a, b, c)].in_file.count('high_qual') 
		+ tuple_dict[(a, b, c)].in_file.count('low_qual')) 
		for (a, b, c) in tuple_dict])
	l = [x for y in l for x in y]  #merge sublists
	max_count = max([l.count(x) for x in l])
	plur = list(set([x for x in l if l.count(x) == max_count]))
	return ([] if len(plur) == len(set(l)) else plur)
		
		
def min_chromosome(chroms):
	"""Find name of chromosome which would appear first in ordered file.
	
	Args: 
		chroms (list of str): Chromosome names ['chrN', 'chrP', ...].
	
	Returns: 
		str: Name of chromosome which would appear first in ordered file.
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
	"""Update counts map given a VCF record. Assumes high-quality variant.
	
	Args:
		counts (dict): Dict of category:count.
		reader (VCF reader): VCF reader for file.
		record (VCF reader record): Specific VCF record.
	"""		
	counts['High-quality variants'] += 1
	counts['Heterozygous'] += record.num_het
	counts['Homozygous alternate'] += record.num_hom_alt
	counts['Homozygous reference'] += record.num_hom_ref
	counts['Missing'] += record.num_unknown
	counts['Indel'] += record.is_indel
	is_snp = (
		1 if (record.is_snp and record.genotype(reader.samples[0]).is_variant) 
		else 0)
	counts['SNP'] += is_snp
	if is_snp:
		counts[
			'%s>%s' % (record.alleles[0], record.alleles[int(
			 max(record.genotype(reader.samples[0]).gt_alleles))])
			] += 1
	

def print_dict_list(
	out_f, table_title, dict_list, cats, header_list, table_message):
	"""Create HTML file and print a list of category:value dictionaries as an
	HTML table.
	
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
		'%dpx; word-wrap:break-word"><tr><th></th>' % (175 * 
		(len(dict_list) + 1)))
	for x in header_list:
		out_f.write('<th>%s</th>' % x)
	out_f.write('</tr>')
	for cat in cats:
		out_f.write('<tr><th>%s</th>' % cat)
		for d in dict_list:
			out_f.write('<td>%s</td>'% (d[cat]))
		out_f.write('</tr>')	
	out_f.write('</table>')
		
		
def normalize_counts(categories, counts):
	"""Add normalized values to counts dict. Used for displaying
	HTML table.
	"""
	if len(counts) == 0:
		return []
	min_index = 0
	flag_dist = .05  #normalized table entries at least flag_dist
					 #distance from 1 are flagged (colored).
	counts_norm = [dict(z) for z in counts]
	for cat in categories:
		min_val = counts[min_index][cat]
		for f in xrange(len(counts)):
			if (
				all(z[cat] == 0 and (not isinstance(z[cat], bool)) 
				for z in counts)):
				counts_norm[f][cat] = str(counts[f][cat]) + ' (%.2f)' % 1
			else:
				if (
					(not isinstance(counts[f][cat], int)) or (isinstance(
					counts[f][cat], bool)) or (min_val == 0)):
					continue
				counts_norm[f][cat] = (
					(counts[f][cat] / float(min_val) * 
					(float(counts[min_index]['High-quality variants']) / 
					counts[f]['High-quality variants'])))
				if abs(1.0 - counts_norm[f][cat]) < flag_dist:
					counts_norm[f][cat] = (
						'%d (%.2f)' % (counts[f][cat], counts_norm[f][cat]))
				else:
					counts_norm[f][cat] = (
						'%d (<font color = "blue">%.2f</font>)' 
						% (counts[f][cat], counts_norm[f][cat]))
	return counts_norm


def output_qc_plot(
	qc_categories, qc_param_label, qc_counts, out_filename, intervals):
	"""Create QC plot using pyplot. Save plot to out_filename.
	"""
	plt.clf()
	for categ in qc_categories:
		if categ != 'all':
			plt.plot(
				intervals, [(1.0 - (float(qc_counts[categ][n])  /
				qc_counts[categ][0])) for n in intervals], label = categ)
	plt.axis([0, max(intervals), 0, 1.05])
	plt.xlabel(qc_param_label)
	plt.ylabel('Fraction removed')
	plt.title('%s as a QC parameter' % qc_param_label)
	ax = plt.subplot(111)
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * .75, box.height])
	ax.legend(
		[z for z in qc_categories if z != 'all'], loc = 'center left', 
		bbox_to_anchor = (1, .87))
	plt.savefig(out_filename)
	
if __name__ == '__main__':
	main()
