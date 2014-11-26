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
		'--out_interesting_variants', type = str, 
		default = 'interesting_variants.txt')	#TODO: Implement this.
	parser.add_argument('vcf_file_names', nargs = '+', type = str)
	args = parser.parse_args()

	in_files = []	# input VCF files
	counts = [] 	# category:count map for each file
	readers = []	# PyVCF reader for each file
	curr_recs = [] 	# current record for each file
	done_readers = []	# indices of readers that have reached end of file
	#statistics, in the order in which they will be displayed
	categories=[
		'Filename', 'Multi_sample', 'Records_Total', 'Records_Skipped', 
		'Hom_Ref', 'SNP', 'Het', 'Hom_Alt', 'Indel', 'Missing', 'A>C', 'A>G', 
		'A>T', 'C>A', 'C>G', 'C>T', 'G>A', 'G>C', 'G>T', 'T>A', 'T>C', 'T>G']	
			
	init_readers(
	args.vcf_file_names, counts, in_files, curr_recs, categories, readers)
		
	with open(args.out_all_variants, 'w') as out_all_variants:
		out_all_variants.write(
			'Chr\tPos\t' + '\t'.join([z['Filename'] for z in counts]) + '\n')		
		compare_variants(
			readers, curr_recs, done_readers, out_all_variants, counts)

	for f in in_files:
		f.close()	
	
	print_dict_list(
		args.out_stats, 'VCF Comparison', 'Counts', counts, categories, 
		'Filename')


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
		counts[-1]['Multi_sample'] = len(reader.samples) > 1
		for cat in [z for z in categories if z not in 
		['Filename', 'Multi_sample']]:
			counts[-1][cat] = 0 if ((not counts[-1]['Multi_sample']) 
			or ('Records' in cat)) else '--'
		readers.append(reader)
		record = reader.next()
		curr_recs.append(record)
		update_counts(counts[-1], reader, record)

	
def compare_variants(
	readers, curr_recs, done_readers, out_all_variants, counts):
	"""Process variants present in VCF files: Iterate through all files 
	simultaneously. Write per-file variant info to all_variants output file 
	and update overall statistics (counts). 
	
	Args:
		readers (list of vcf.Reader): VCF readers.
		curr_recs (list of vcf.model._Record): Record that each file is on.
		done_readers (list of int): 
			Indices of readers that have reached end of file.
		out_all_variants (text file): 
			Open output file to write all variants to.
		counts (list of dicts): Category:count dict per file.
	"""
	
	min_chrom = min_chromosome(
		list(curr_recs[z].CHROM for z in xrange(len(curr_recs)) 
		if not z in done_readers))
	print 'Processing ' + min_chrom + '...'
	while len(done_readers) < len(readers):
		new_min_chrom = min_chromosome(
			list(curr_recs[z].CHROM for z in xrange(len(curr_recs)) 
			if not z in done_readers))
		if new_min_chrom != min_chrom:
			print 'Processing ' + new_min_chrom + '...'
			min_chrom = new_min_chrom
		min_pos = min(
			curr_recs[w].POS for w in xrange(len(curr_recs)) 
			if (curr_recs[w].CHROM == min_chrom and not w in done_readers))
		out_all_variants.write(min_chrom + ' ' + str(min_pos) + '\t')
		for i, reader in enumerate(readers):
			if (
				curr_recs[i].CHROM == min_chrom
				and curr_recs[i].POS == min_pos and not i in done_readers):
				out_all_variants.write('Y\t')
				try:
					record = reader.next()
					update_counts(counts[i], reader, record)
					curr_recs[i] = record
				except StopIteration:
					done_readers.append(i)
			else:
				out_all_variants.write('N\t')		
		out_all_variants.write('\n')


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
	counts['Records_Total'] += 1	
	if counts['Multi_sample']:
		counts['Records_Skipped'] += 1
		return		
	counts['Het'] += record.num_het
	counts['Hom_Alt'] += record.num_hom_alt
	counts['Hom_Ref'] += record.num_hom_ref
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
	out_f_name, page_title, table_title, dict_list, cats, header_cat):
	"""Create HTML file and print a list of category:value dictionaries as an
	HTML table.
	
	Args:
		out_f_name (str): File path to write to.
		page_title (str): HTML page title.
		table_title (str): HTML table title.
		dict_list (list of dict): Category:value dicts.
		cats (list of string): Ordered categories for stats output.
		header_cat (str): Category to group other categories by.
	"""
	with open(out_f_name, 'w') as out_f:
		out_f.write(
			'<!DOCTYPE html><html><head><title>' + page_title 
			+ '</title><body><h1>' + page_title + '</h1>')
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
		out_f.write('</body></html>')


if __name__ == '__main__':
	main()