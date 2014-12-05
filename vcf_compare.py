"""Compares VCF files. Produces statistics and variant comparison files."""

import os
import sys
import argparse
import vcf
import itertools

def main():
	parser = argparse.ArgumentParser()	
	parser.add_argument('--out_stats', type = str, default = 'vcf_stats.html')
	parser.add_argument(
		'--out_all_variants', type = str, default = 'all_variants.txt')
	parser.add_argument(
		'--out_interesting_variants', type = str, 
		default = 'interesting_variants.txt')
	#min quality, depth for inclusion in interesting variants
	parser.add_argument('--min_qual', type = int, default = 20) 
	parser.add_argument('--min_depth', type = int, default = 30)
	parser.add_argument('vcf_file_names', nargs = '+', type = str)

	args = parser.parse_args()

	in_files = []	# input VCF files
	counts = [] 	# category:count map for each file
	readers = []	# PyVCF reader for each file
	curr_recs = [] 	# current record for each file
	#statistics, in the order in which they will be displayed
	categories=[
		'Filename', 'Total Variants', 'Homozygous Reference', 'SNP', 
		'Heterozygous', 'Homozygous Alternate', 'Indel', 'Missing', 
		'A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T', 'G>A', 'G>C', 'G>T', 
		'T>A', 'T>C', 'T>G', 'Multi-sample', 'Records Skipped']
			
	init_readers(
		args.vcf_file_names, counts, in_files, curr_recs, categories, readers)

	header = (
		'Chr\tPos\tRef\tAlt\t' + '\t'.join([z['Filename'] for z in counts]) 
		+ '\t%Y\tAve_GQ\tAve_DP' + '\n')
	with open(args.out_all_variants, 'w') as out_all_variants, \
	open(args.out_interesting_variants, 'w') as out_interesting_variants:	
		out_all_variants.write(header)		
		out_interesting_variants.write(header)
		compare_variants(
			readers, curr_recs, out_all_variants, out_interesting_variants, 
			counts, args.min_qual, args.min_depth)

	for f in in_files:
		f.close()	
	
	#make list of normalized counts dicts
	min_index = (
		[z['Total Variants'] for z in counts].index(min([y['Total Variants'] 
		for y in counts if not y['Multi-sample']])))
	counts_norm = [dict(z) for z in counts]
	for cat in categories:
		min_val = counts[min_index][cat]
		for f in xrange(len(counts)):	
			if isinstance(counts[f][cat], int) and min_val != 0:
				counts_norm[f][cat] = (
					str(counts[f][cat]) + ' (%.2f' % (counts[f][cat] / 
					float(min_val)) + ')')
	
	with open(args.out_stats, 'w') as out_f:
		page_title = 'VCF Comparison'
		out_f.write(
			'<!DOCTYPE html><html><head><title>' + page_title 
			+ '</title><body><h1>' + page_title + '</h1>')
# 		print_dict_list(
# 			out_f, 'Counts', counts, categories, 'Filename')
		print_dict_list(
			out_f, 'Counts (with normalization)', counts_norm, categories,
			'Filename')
		out_f.write('</body></html>')

def init_readers(f_names, counts, in_files, curr_recs, categories, readers):
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
				'Records Skipped') or (cat == 'Total Variants')) else '--')
		readers.append(reader)
		record = reader.next()
		curr_recs.append(record)
		update_counts(counts[-1], reader, record)

	
def compare_variants(
	readers, curr_recs, out_all_variants, out_interesting_variants, counts, 
	min_qual, min_depth):
	"""Process variants present in VCF files: Iterate through all files 
	simultaneously. Write per-file variant info to all_variants output file 
	and update overall statistics (counts). 
	
	Args:
		readers (list of vcf.Reader): VCF readers.
		curr_recs (list of vcf.model._Record): Record that each file is on.
		out_all_variants (text file): 
			Open output file to write all variants to.
		out_interesting_variants (text file): 
			Open output file to write disagreements above min_qual and 
			min_depth to.
		counts (list of dicts): Category:count dict per file.
		min_qual (int): Min GQ score for inclusion in interesting_variants.
		min_depth (int): Min DP score for inclusion in interesting_variants.
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
			if (
				curr_recs[i].CHROM == min_chrom
				and curr_recs[i].POS == min_pos and not i in done_readers):
				curr_variant_in_file.append(True)
				curr_calls[i] = curr_recs[i].genotype(reader.samples[0])
				try:
					record = reader.next()
					update_counts(counts[i], reader, record)
					curr_recs[i] = record
				except StopIteration:
					done_readers.append(i)
			else:
				curr_variant_in_file.append(False)

		output_variant(
			curr_variant_in_file, prev_recs, curr_calls, out_all_variants,
			out_interesting_variants, min_qual, min_depth)

class RefAlt:
	def __init__(self, num_records):
		self.total_gq = 0
		self.count_gq = 0
		self.total_dp = 0
		self.count_dp = 0
		self.in_file = [False] * num_records
		
def output_variant(
	variant_in_file, records, calls, out_all_variants, 
	out_interesting_variants, min_qual, min_depth):
	"""Output a single variant to out_all_variants and out_interesting_variants.
	
		Args:
			variant_in_file (list of bool): True iff file has variant.
			records (list of vcf.model._Record): Record of variant, per file.
			calls (list of vcf.model._Call): Variant call, per file.
			out_all_variants (text file): 
				Open output file to write all variants to.
			out_interesting_variants (text file): 
				Open output file to write disagreements above min_qual and 
				min_depth to.
			min_qual (int): Min GQ score for inclusion in interesting_variants.
			min_depth (int): Min DP score for inclusion in interesting_variants.
	"""
	#high-quality variants that were present in some files and not others
	
	#print records[i] 	Record(CHROM=chr1, POS=900001, REF=G, ALT=[A])
	#print calls[i]		Call(sample=1463-02-p1-2500rapid, CallData(GT=0/1, AD=[3, 3], DP=6, GQ=82.25, PL=[82, 0, 94]))

	high_qual_var_in_file = False
	variant_in_all_files = True
	multiple_high_qual_vars = False
	
	dict_var_stats = {} 	#(ref, alt) : RefAlt obj
	for i in range(len(records)):
		if variant_in_file[i]:
			ref = str(records[i].REF)
			alt = str(records[i].ALT[0])
			new_var = False
			if not (ref, alt) in dict_var_stats:
				new_var = True
				dict_var_stats[(ref, alt)] = RefAlt(len(records))
			dict_var_stats[(ref, alt)].in_file[i] = True
			gq = calls[i]['GQ']
			dp = calls[i]['DP']
			if gq is not None:
				dict_var_stats[(ref, alt)].total_gq += gq
				dict_var_stats[(ref, alt)].count_gq += 1
			if dp is not None:
				dict_var_stats[(ref, alt)].total_dp += dp
				dict_var_stats[(ref, alt)].count_dp += 1
			if (
				gq is not None and dp is not None and gq >= min_qual and 
				dp >= min_depth):
				if high_qual_var_in_file and new_var:
					multiple_high_qual_vars = True
				high_qual_var_in_file = True				
		else:
			variant_in_all_files = False
	
	chrom = records[variant_in_file.index(True)].CHROM
	pos = records[variant_in_file.index(True)].POS

	for x in dict_var_stats:
		ref = x[0]
		alt = x[1]
		ave_gq = float(dict_var_stats[x].total_gq) / dict_var_stats[x].count_gq
		ave_dp = float(dict_var_stats[x].total_dp) / dict_var_stats[x].count_dp
		percent_y = (
			100.0 * dict_var_stats[(ref, alt)].in_file.count(True) 
			/ len(records))
		variant_in_file_str = (
			'\t'.join('Y' if z else 'N' for z in 
			dict_var_stats[(ref, alt)].in_file))
		
		out_str = (
			'%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\n' % (chrom, 
			str(pos), ref, alt, variant_in_file_str, percent_y, 
			ave_gq, ave_dp))
		out_all_variants.write(out_str)	
		
		if (
			(high_qual_var_in_file and not variant_in_all_files) or 
			multiple_high_qual_vars):
			out_interesting_variants.write(out_str)


def min_chromosome(chroms):
	"""Find name of chromosome which would appear first in ordered file.
	
	Args: 
		chroms (list of str): Chromosome names ['chrN', 'chrP', ...].
	
	Returns: 
		str: Name of chromosome which would appear first in ordered file.
	"""
	chroms_processed = list(chroms)
	for i in range(len(chroms_processed)):
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
	counts['Total Variants'] += 1	
	if counts['Multi-sample']:
		counts['Records Skipped'] += 1
		return		
	counts['Heterozygous'] += record.num_het
	counts['Homozygous Alternate'] += record.num_hom_alt
	counts['Homozygous Reference'] += record.num_hom_ref
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
	out_f, table_title, dict_list, cats, header_cat):
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
		'<h2>' + table_title + '</h2><table border="1"><tr><th></th>')
	for d in dict_list:
		out_f.write('<th>' + d[header_cat] + '</th>')	
	out_f.write('</tr>')
	for cat in [z for z in cats if z != header_cat]:
		out_f.write('<tr><th>' + cat + '</th>')
		for d in dict_list:
			out_f.write('<td>' + str(d[cat]) + '</td>')
		out_f.write('</tr>')	
	out_f.write('</table>')


if __name__ == '__main__':
	main()
