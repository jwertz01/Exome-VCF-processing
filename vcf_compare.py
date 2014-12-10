"""Compares VCF files. Produces statistics and variant comparison files."""

import os
import sys
import argparse
import vcf

def main():
	parser = argparse.ArgumentParser()	
	parser.add_argument('--out_stats', type = str, default = 'vcf_stats.html')
	parser.add_argument(
		'--out_all_variants', type = str, default = 'all_variants.txt')
	parser.add_argument(
		'--out_conflicting_variants', type = str, 
		default = 'conflicting_variants.txt')
	#min quality, depth for inclusion in variant output files
	parser.add_argument('--min_qd', type = int, default = 3) 
	parser.add_argument('--min_dp', type = int, default = 20)
	parser.add_argument('vcf_file_names', nargs = '+', type = str)

	args = parser.parse_args()

	in_files = []	# input VCF files
	counts = [] 	# category:count map for each file
	readers = []	# PyVCF reader for each file
	curr_recs = [] 	# current record for each file
	#statistics, in the order in which they will be displayed
	categories=[
		'Filename', 'Total variants', 'Variants that pass quality filter',
		'Homozygous reference', 'SNP', 'Heterozygous', 'Homozygous alternate',
		'Indel', 'Missing', 'A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T', 'G>A', 
		'G>C', 'G>T', 'T>A', 'T>C', 'T>G', 'Multi-sample', 'Records skipped',
		'Variants absent from this file present in other files',
		'Variants present that conflict with other files']
			
	init_readers(
		args.vcf_file_names, counts, in_files, curr_recs, categories, readers,
		args.min_qd, args.min_dp)

	header = (
		'Chr\tPos\tRef\tAlt\t' + '\t'.join([z['Filename'] for z in counts]) 
		+ '\t%Present\tAve_QD\tAve_DP' + '\n')
		
	with open(args.out_all_variants, 'w') as out_all_variants, \
	open(args.out_conflicting_variants, 'w') as out_conflicting_variants:	
		out_all_variants.write(header)		
		out_conflicting_variants.write(header)
		compare_variants(
			readers, curr_recs, out_all_variants, out_conflicting_variants, 
			counts, args.min_qd, args.min_dp)

	for f in in_files:
		f.close()	
	
	counts_norm = normalize_counts(categories, counts)
	with open(args.out_stats, 'w') as out_f:
		message = (
			'Note: Normalized with respect to the number of variants that ' +
			'pass the quality filter, and to the value in the leftmost '+
			'numeric column. Normalized values less than 0.95 or ' +
			'greater than 1.05 are shown in blue. Counts shown are after ' +
			'quality filtering.')
		print_dict_list(
			out_f, 'VCF Comparison', 'Counts (with normalization)', 
			counts_norm, categories, 'Filename', message)


def init_readers(
	f_names, counts, in_files, curr_recs, categories, readers, min_qd, min_dp):
	"""
	Initialize VCF readers and stats dict.
	
	Args:
		f_names (list of str): VCF file names.
		counts (list of dict): Category:count dict per file.
		in_files (list of file): All VCF files.
		curr_recs (list of vcf.model._Record): Record that each file is on.
		categories (list of str): Ordered categories for stats output.
		readers (list of vcf.Reader): VCF readers.
	"""
	for f_name in f_names:
		in_files.append(open(f_name))
		counts.append({})
		f_base_name = os.path.basename(f_name)
		if '.' in f_base_name:
			f_base_name = f_base_name[:f_base_name.index('.')]		
		counts[-1]['Filename'] = f_base_name
		reader = vcf.Reader(in_files[-1])
		counts[-1]['Multi-sample'] = len(reader.samples) > 1
		for cat in [z for z in categories if z not in 
		['Filename', 'Multi-sample']]:
			counts[-1][cat] = (
				0 if ((not counts[-1]['Multi-sample']) or (cat == 
				'Records skipped') or (cat == 'Total variants')) else '--')
		readers.append(reader)
		is_high_qual = False
		record = reader.next()
		curr_recs.append(record)
		

def compare_variants(
	readers, curr_recs, out_all_variants, out_conflicting_variants, counts, 
	min_qd, min_dp):
	"""Process variants present in VCF files: Iterate through all files 
	simultaneously. Write per-file variant info to all_variants output file 
	and update overall statistics (counts). 
	
	Args:
		readers (list of vcf.Reader): VCF readers.
		curr_recs (list of vcf.model._Record): Record that each file is on.
		out_all_variants (text file): 
			Open output file to write all variants above min_qd and min_dp
		 	to.
		out_conflicting_variants (text file): 
			Open output file to write disagreements above min_qd and 
			min_dp to.
		counts (list of dicts): Category:count dict per file.
		min_qd (int): Min QD score for inclusion in variant output files.
		min_dp (int): Min DP score for inclusion in variant output files.
	"""
	done_readers = []
	min_chrom = min_chromosome(
		list(curr_recs[z].CHROM for z in xrange(len(curr_recs)) 
		if not z in done_readers))
	print 'Processing ' + min_chrom + '...'
	while len(done_readers) < len(readers):
		curr_variant_in_file = []
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
			call = curr_recs[i].genotype(reader.samples[0])
			qd = curr_recs[i].INFO['QD']
			dp = call['DP']
			if (
				curr_recs[i].CHROM == min_chrom and 
				curr_recs[i].POS == min_pos and not i in done_readers):
				counts[i]['Total variants'] += 1
				is_high_qual = (
					dp is not None and qd is not None and qd >= min_qd and 
					dp >= min_dp)
				if is_high_qual:
					curr_variant_in_file.append('high_qual')
					curr_calls[i] = call
					update_counts(counts[i], reader, curr_recs[i])
				else:
					curr_variant_in_file.append('low_qual')
				try:
					record = reader.next()
					curr_recs[i] = record
				except StopIteration:
					done_readers.append(i)
			else:
				curr_variant_in_file.append('absent')
		output_variant(
			curr_variant_in_file, prev_recs, curr_calls, out_all_variants,
			out_conflicting_variants, min_qd, min_dp, counts)


#used with output_variant function		
class RefAlt:
	def __init__(self, num_records):
		self.total_qd = 0
		self.count_qd = 0
		self.total_dp = 0
		self.count_dp = 0
		self.in_file = ['absent'] * num_records
		
def output_variant(
	variant_in_file, records, calls, out_all_variants, 
	out_conflicting_variants, min_qd, min_dp, counts):
	"""Output a single variant to out_all_variants and out_conflicting_variants.
	
		Args:
			variant_in_file (list of bool): True iff file has variant.
			records (list of vcf.model._Record): Record of variant, per file.
			calls (list of vcf.model._Call): Variant call, per file.
			out_all_variants (text file): 
				Open output file to write all variants above min_qd and 
				min_dp to.
			out_conflicting_variants (text file): 
				Open output file to write disagreements above min_qd and 
				min_dp to.
			min_qd (int): Min QD score for inclusion in variant output files.
			min_dp (int): Min DP score for inclusion in variant output files.
			counts (list of dicts): Category:count dict per file.
	"""
	if all(z in ['low_qual', 'absent'] for z in variant_in_file):
		return
	variant_in_all_files = True
		
	chrom = records[variant_in_file.index('high_qual')].CHROM
	pos = records[variant_in_file.index('high_qual')].POS
	
	#find stats for each different (ref, alt) combination at given position
	#on chromosome
	dict_var_stats = {} 	#(ref, alt) : RefAlt obj
	for i in xrange(len(records)):
		if variant_in_file[i] == 'absent':
			continue
		ref = str(records[i].REF)
		alt = str(records[i].ALT[0])
		if not (ref, alt) in dict_var_stats:
			dict_var_stats[(ref, alt)] = RefAlt(len(records))
		dict_var_stats[(ref, alt)].in_file[i] = variant_in_file[i]
		if variant_in_file[i] == 'high_qual':
			qd = records[i].INFO['QD']
			dp = calls[i]['DP']
			if qd is not None:
				dict_var_stats[(ref, alt)].total_qd += qd
				dict_var_stats[(ref, alt)].count_qd += 1
			if dp is not None:
				dict_var_stats[(ref, alt)].total_dp += dp
				dict_var_stats[(ref, alt)].count_dp += 1				
			
	#output each (ref, alt) combination at given position on chromosome
	for x in dict_var_stats:
		if all (
			z != 'high_qual' for z in dict_var_stats[x].in_file):
			continue
		ref = x[0]
		alt = x[1]
		ave_qd = float(dict_var_stats[x].total_qd) / dict_var_stats[x].count_qd
		ave_dp = float(dict_var_stats[x].total_dp) / dict_var_stats[x].count_dp
		percent_present = (
			100.0 * (dict_var_stats[x].in_file.count('low_qual') + 
			dict_var_stats[x].in_file.count('high_qual')) / len(records))
		variant_in_file_str = '\t'.join(dict_var_stats[x].in_file)
		for i, in_file in enumerate(dict_var_stats[x].in_file):
			if not counts[i]['Multi-sample']:
				if in_file == 'absent':
					#same chrom, pos but different ref/ alt counts as absent.
					counts[i]['Variants absent from this file ' +
					'present in other files'] += 1
					counts[i]['Variants present that conflict ' +
					'with other files'] += 1
		out_str = (
			'%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\n' % (chrom, 
			str(pos), ref, alt, variant_in_file_str, percent_present, 
			ave_qd, ave_dp))
		out_all_variants.write(out_str)	
		if percent_present not in [0, 100]:
			out_conflicting_variants.write(out_str)

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
	"""Update counts map given a VCF record. 
	
	Args:
		counts (dict): Dict of category:count.
		reader (VCF reader): VCF reader for file.
		record (VCF reader record): Specific VCF record.
	"""	
	if counts['Multi-sample']:
		counts['Records skipped'] += 1
		return		
	counts['Variants that pass quality filter'] += 1
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
	out_f, page_title, table_title, dict_list, cats, header_cat, message):
	"""Create HTML file and print a list of category:value dictionaries as an
	HTML table.
	
	Args:
		out_f (str): File to write to.
		table_title (str): HTML table title.
		dict_list (list of dict): Category:value dicts.
		cats (list of string): Ordered categories for stats output.
		header_cat (str): Category to group other categories by.
	"""
	out_f.write(
		'<!DOCTYPE html><html><head><title>' + page_title 
		+ '</title></head><body><h1>' + page_title + '</h1>')
	out_f.write(
		'<h2>' + table_title + '</h2><table border="1"><tr><th></th>')
	for d in dict_list:
		out_f.write(
			'<th>' + d[header_cat] + '</th>')
	out_f.write('</tr>')
	for cat in [z for z in cats if z != header_cat]:
		out_f.write(
			'<tr><th style="width:225px">' + cat + '</th>')
		for d in dict_list:
			out_f.write('<td>' + str(d[cat]) + '</td>')
		out_f.write('</tr>')	
	out_f.write('</table>')
	out_f.write('<p>' + message + '</p>')
	out_f.write('</body></html>')
		
		
def normalize_counts(categories, counts):
	min_index = (
		min(i for i in xrange(len(counts)) if not counts[i]['Multi-sample']))
	counts_norm = [dict(z) for z in counts]
	for cat in categories:
		if (
			all((isinstance(z[cat], int) and (not isinstance(z[cat], bool))
			and z[cat] == 0) for z in counts)):
			for f in xrange(len(counts)):
				counts_norm[f][cat] = str(counts[f][cat]) + ' (%.2f)' % 1
		else:
			min_val = counts[min_index][cat]
			for f in xrange(len(counts)):	
				if (
					(not isinstance(counts[f][cat], int)) or (isinstance(
					counts[f][cat], bool)) or (min_val == 0) or
					(counts[f]['Multi-sample'])):
					continue
				counts_norm[f][cat] = (
					(counts[f][cat] / float(min_val) * 
					(float(counts[min_index]['Variants that pass ' +
					'quality filter']) / 
					counts[f]['Variants that pass quality filter'])))
				if (
					counts_norm[f][cat] < 1.05 and counts_norm[f][cat] 
					> 0.95):
					counts_norm[f][cat] = (
						str(counts[f][cat]) + ' (%.2f)' % counts_norm[f][cat])
				else:
					counts_norm[f][cat] = (
						str(counts[f][cat]) + 
						' (<font color = "blue">%.2f</font>)' 
						% counts_norm[f][cat])	
	return counts_norm

if __name__ == '__main__':
	main()
