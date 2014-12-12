"""Compares VCF files. Produces statistics and variant comparison files."""

import os
import sys
import argparse
import vcf
import matplotlib.pyplot as plt

def main():
	parser = argparse.ArgumentParser()	
	parser.add_argument('vcf_file_names', nargs = '+', type = str)
	parser.add_argument('--out_stats', type = str, default = 'vcf_stats.html')
	parser.add_argument(
		'--depth_plot', type = str, default = 'dp_plot.png')
	parser.add_argument(
		'--qd_plot', type = str, default = 'qd_plot.png')
	parser.add_argument(
		'--qual_plot', type = str, default = 'qual_plot.png')
	parser.add_argument('--create_qc_plots', type = str)	
	# containing all high-quality variants
	parser.add_argument(	
		'--out_all_variants', type = str, default = 'all_variants.txt')
		
	# containing variants such that at least one file has high-quality 
	# variant and another file does not have variant or has different 
	# high-quality ref/alt combination at chrom and pos.
	parser.add_argument(
		'--out_conflicting_variants', type = str, default = 
		'conflicting_variants.txt')
	
	# containing variants such that at least one file has high-quality
	# variant and another file has different high-quality ref/alt combination 
	# at chrom and pos.
	parser.add_argument(
		'--out_interesting_variants', type = str, 
		default = 'interesting_variants.txt')
		
	# min quality by depth, read depth for inclusion 
	# in above variant output files
	parser.add_argument('--min_qd', type = int, default = 7) #change to 5-10
	parser.add_argument('--min_dp', type = int, default = 30)

	args = parser.parse_args()

	in_files = []	# input VCF files
	multi_sample_files = []	#names of multi-sample VCF files
	counts = [] 	# category:count map for each file
	readers = []	# PyVCF reader for each file
	curr_recs = [] 	# current record for each file
	#qc plot x-axis intervals
	qc_intervals = {
		'dp': xrange(0, 55, 5), 'qd': xrange(0, 20, 2), 
		'qual': xrange(0, 260, 26)}
	qc_categories = [
		'all', 'het_all_agree', 'het_disagree', 'hom_all_agree', 
		'hom_disagree']
	qc_counts = {}
	
	#statistics, in the order in which they will be displayed
	#in html table
	categories=[
		'Filename', 'Total variants', 'High-quality variants',
		'Homozygous reference', 'SNP', 'Heterozygous', 'Homozygous alternate',
		'Indel', 'Missing', 'A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T', 'G>A', 
		'G>C', 'G>T', 'T>A', 'T>C', 'T>G',
		'Variants absent from this file present in other files',
		'High-quality variants present that conflict with other files']
			
	#initialize readers, counts
	for f_name in args.vcf_file_names:
		in_f = open(f_name)
		reader = vcf.Reader(in_f)
		f_base_name = os.path.basename(f_name)		
		if len(reader.samples) > 1:
			multi_sample_files.append(f_base_name)
			in_f.close()
		else:
			in_files.append(in_f)
			counts.append({})
			for cat in categories:
				counts[-1][cat] = 0
			counts[-1]['Filename'] = f_base_name
			readers.append(reader)
			curr_recs.append(reader.next())
	
	for i in qc_intervals:
		qc_counts[i] = {}			
		for j in qc_categories:
			qc_counts[i][j] = {}
			for k in qc_intervals[i]:
				qc_counts[i][j][k] = 0	
			
	#iterate over vcf files, write variant comparison output files, 
	#fill in counts	
	header = (
		'Chr\tPos\tRef\tAlt\t' + '\t'.join([z['Filename'] for z in counts]) 
		+ '\t%Present\tAve_QD\tAve_DP' + '\n')
	with open(args.out_all_variants, 'w') as out_all_variants, \
	open(args.out_conflicting_variants, 'w') as out_conflicting_variants, \
	open(args.out_interesting_variants, 'w') as out_interesting_variants:	
		out_all_variants.write(header)		
		out_conflicting_variants.write(header)
		out_interesting_variants.write(header)
		compare_variants(
			readers, curr_recs, out_all_variants, out_conflicting_variants, 
			out_interesting_variants, counts, args.min_qd, args.min_dp,
			qc_counts)
	for f in in_files:
		f.close()	
	
	#output qc plots
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
 	counts_norm = normalize_counts(categories, counts)
	with open(args.out_stats, 'w') as out_f:
		page_title = 'VCF Comparison'
		out_f.write(
			'<!DOCTYPE html><html><head><title>' + page_title 
			+ '</title></head><body><h1>' + page_title + '</h1>')
		table_message = (
			'Note: Normalized with respect to the number of ' +
			'high-quality variants, and to the value in the leftmost ' +
			'numeric column. Normalized values less than 0.95 or ' +
			'greater than 1.05 are shown in blue. Counts shown are after ' +
			'quality filtering.')
		print_dict_list(
			out_f, 'Counts (with normalization)', 
			counts_norm, categories, 'Filename', table_message)
		if len(multi_sample_files) > 0:
			out_f.write(
				'<p><b>Multi-sample files excluded from ' +
				'analysis: </b>' + ', '.join(multi_sample_files) + '</p>')
		out_f.write('</body></html>')
		
		
def compare_variants(
	readers, curr_recs, out_all_variants, out_conflicting_variants, 
	out_interesting_variants, counts, min_qd, min_dp, qc_counts):
	"""Iterate through all VCF files simultaneously. Write per-file 
	variant info to output files and update overall statistics (counts). 
	"""
	if min(len(readers), len(curr_recs), len(counts)) == 0:
		return
	done_readers = [] #indices of readers for which eof has been reached
			
	min_chrom = min_chromosome(
		list(curr_recs[z].CHROM for z in xrange(len(curr_recs)) 
		if not z in done_readers))
	print 'Processing ' + min_chrom + '...'
	while len(done_readers) < len(readers):
		curr_variant_in_file = [] # bool list. True iff file has variant.
		new_min_chrom = min_chromosome(
			list(curr_recs[z].CHROM for z in xrange(len(curr_recs)) 
			if not z in done_readers))
		if new_min_chrom != min_chrom:
			print 'Processing ' + new_min_chrom + '...'
			min_chrom = new_min_chrom
		min_pos = min(
			curr_recs[w].POS for w in xrange(len(curr_recs)) 
			if (curr_recs[w].CHROM == min_chrom and not w in done_readers))
		curr_calls = [None] * len(readers)
		prev_recs = list(curr_recs)
		for i, reader in enumerate(readers):
			record = curr_recs[i]
			call = record.genotype(reader.samples[0])
			qd = record.INFO['QD']
			dp = call['DP']
			if (
				record.CHROM == min_chrom and 
				record.POS == min_pos and not i in done_readers):
				counts[i]['Total variants'] += 1		
				is_high_qual = (
					dp is not None and qd is not None and qd >= min_qd and 
					dp >= min_dp)
				if is_high_qual:
					curr_variant_in_file.append('high_qual')
					update_counts(counts[i], reader, curr_recs[i])
				else:
					curr_variant_in_file.append('low_qual')
				curr_calls[i] = call
				try:
					curr_recs[i] = reader.next()
				except StopIteration:
					done_readers.append(i)
			else:
				curr_variant_in_file.append('absent')
				
		#write single variant (at min_chrom and min_pos) to output files,
		#update counts
		process_variant(
			curr_variant_in_file, prev_recs, curr_calls, out_all_variants,
			out_conflicting_variants, out_interesting_variants, min_qd, min_dp,
			counts, qc_counts)

	
#used with process_variant function		
class RefAlt:
	def __init__(self, num_records):
		self.total_qd = 0
		self.count_qd = 0
		self.total_dp = 0
		self.count_dp = 0
		self.in_file = ['absent'] * num_records
		
		
def process_variant(
	variant_in_file, records, calls, out_all_variants, 
	out_conflicting_variants, out_interesting_variants, min_qd, min_dp, 
	counts, qc_counts):
	"""Output a single variant to output files, and update counts.
	"""
	high_qual_present = (
		not all(z in ['low_qual', 'absent'] for z in variant_in_file))
	
	#find each (ref, alt) combination at given position
	#on chromosome
	dict_var_stats = {} 	#(ref, alt) : RefAlt obj

	for i in xrange(len(records)):
		if variant_in_file[i] == 'absent':
			continue		
		call = calls[i]
		record = records[i]
		ref = str(record.REF)
		alt = str(record.ALT[0])
		try:
			dp = call['DP']
		except TypeError:
			dp = None 
		qd = record.INFO['QD']
		
		#update dict_var_stats
		if high_qual_present:			
			if not (ref, alt) in dict_var_stats:
				dict_var_stats[(ref, alt)] = RefAlt(len(records))
			dict_var_stats[(ref, alt)].in_file[i] = variant_in_file[i]									
			if variant_in_file[i] == 'high_qual':
				if qd is not None:
					dict_var_stats[(ref, alt)].total_qd += qd
					dict_var_stats[(ref, alt)].count_qd += 1
				if dp is not None:
					dict_var_stats[(ref, alt)].total_dp += dp
					dict_var_stats[(ref, alt)].count_dp += 1
					
	multiple_vars = (
		len([z for z in dict_var_stats if 'high_qual' in 
		dict_var_stats[z].in_file]) > 1)
		
	#update qc counts
	all_agree = (
		all(v != 'absent' for v in variant_in_file) and (not multiple_vars))
	for i in xrange(len(records)):
		if variant_in_file[i] == 'absent':
			continue
		call = calls[i]
		record = records[i]	
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
		chrom = records[variant_in_file.index('high_qual')].CHROM
		pos = records[variant_in_file.index('high_qual')].POS		
		
		#output each (ref, alt) combination at given position on chromosome
		for x in dict_var_stats:
			 output_ref_alt(x, out_all_variants, out_conflicting_variants, 
			 out_interesting_variants, dict_var_stats, counts, chrom, pos,
			 multiple_vars)
	
	
def update_qc_counts(record, all_agree, qc_val, qc_counts):
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
	ref_alt, out_all_variants, out_conflicting_variants,
	out_interesting_variants, dict_var_stats, counts, chrom, pos,
	multiple_vars):	
	ref_alt_obj = dict_var_stats[ref_alt]
	if all (z != 'high_qual' for z in ref_alt_obj.in_file):
		return
	ref = ref_alt[0]
	alt = ref_alt[1]
	ave_qd = float(ref_alt_obj.total_qd) / ref_alt_obj.count_qd
	ave_dp = float(ref_alt_obj.total_dp) / ref_alt_obj.count_dp
	percent_present = (
		100.0 * (ref_alt_obj.in_file.count('low_qual') + 
		ref_alt_obj.in_file.count('high_qual')) / len(ref_alt_obj.in_file))
	variant_in_file_str = '\t'.join(ref_alt_obj.in_file)
	for i, in_file in enumerate(ref_alt_obj.in_file):
		if in_file == 'absent':
			#same chrom, pos but different ref/ alt counts as absent.
			counts[i]['Variants absent from this file ' +
			'present in other files'] += 1
		if (
			in_file == 'high_qual' and percent_present not in 
			[0, 100]):
			counts[i]['High-quality variants present that ' +
			'conflict with other files'] += 1
	out_str = (
		'%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\n' % (chrom, 
		str(pos), ref, alt, variant_in_file_str, percent_present, 
		ave_qd, ave_dp))
	out_all_variants.write(out_str)	
	if percent_present not in [0, 100]:
		out_conflicting_variants.write(out_str)
		if multiple_vars:
			out_interesting_variants.write(out_str)


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
			record.alleles[0] + '>' + str(record.alleles
			[int(max(record.genotype(reader.samples[0]).gt_alleles))])] += 1
	

def print_dict_list(
	out_f, table_title, dict_list, cats, header_cat, table_message):
	"""Create HTML file and print a list of category:value dictionaries as an
	HTML table.
	
	Args:
		out_f (str): File to write to.
		table_title (str): HTML table title.
		dict_list (list of dict): Category:value dicts.
		cats (list of string): Ordered categories for stats output.
		header_cat (str): Category to group other categories by.
	"""
	if len(dict_list) == 0:
		return
	out_f.write(
		'<h2>' + table_title + '</h2><table border="1"><tr><th></th>')
	for d in dict_list:
		out_f.write('<th>' + d[header_cat] + '</th>')
	out_f.write('</tr>')
	for cat in [z for z in cats if z != header_cat]:
		out_f.write(
			'<tr><th style="width:225px">' + cat + '</th>')
		for d in dict_list:
			out_f.write('<td>' + str(d[cat]) + '</td>')
		out_f.write('</tr>')	
	out_f.write('</table>')
	out_f.write('<p>' + table_message + '</p>')
		
		
def normalize_counts(categories, counts):
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
				if abs(1.0 - counts_norm[f][cat]) < .05:
					counts_norm[f][cat] = (
						str(counts[f][cat]) + ' (%.2f)' % counts_norm[f][cat])
				else:
					counts_norm[f][cat] = (
						str(counts[f][cat]) + 
						' (<font color = "blue">%.2f</font>)' 
						% counts_norm[f][cat])	
	return counts_norm


def output_qc_plot(
	qc_categories, qc_param_label, qc_counts, out_filename, intervals):
	plt.clf()
	for categ in qc_categories:
		if categ != 'all':
			plt.plot(
				intervals, [(1.0 - (float(qc_counts[categ][n])  /
				qc_counts[categ][0])) for n in intervals], label = categ)
	plt.axis([0, max(intervals), 0, 1.05])
	plt.xlabel(qc_param_label)
	plt.ylabel('Fraction removed')
	plt.title(qc_param_label + ' as a QC parameter')
	ax = plt.subplot(111)
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * .75, box.height])
	ax.legend(
		[z for z in qc_categories if z != 'all'], loc = 'center left', 
		bbox_to_anchor = (1, .87))
	plt.savefig(out_filename)
	
	
if __name__ == '__main__':
	main()
