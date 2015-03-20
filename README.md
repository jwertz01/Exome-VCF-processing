# Exome-VCF-processing

##Purpose
Various scripts that compare variants between VCF files. Variant comparison files
and statistics are produced.

vcf_compare.py: Compares any number of VCF files.

vcf_compare_inheritance.py: Additionaly analyzes inheritance.
(Use with VCF files obtained from a parent-offspring trio.)

vcf_compare_exact_same.py: Compares VCFs expected to be
exactly the same as each other. (e.g., the same sample run through
different Galaxy instances, but processed identically otherwise.)

##Requirements
Python 2.x (https://www.python.org/)

PyVCF (https://github.com/jamescasbon/PyVCF)

matplotlib (http://matplotlib.org/)

##Example commands

python vcf_compare.py a.vcf b.vcf --out_stats c.html --out_all_variants
d.txt --out_absent e.txt --out_multiple f.txt --min_qd 5 --min_qual 100
--min_dp 20

python vcf_compare_inheritance.py -p1 a.vcf -p2 b.vcf -c c.vcf
--out_stats d.html --out_all_variants e.txt --min_qd 5 --min_qual 100
--min_dp 20
(Where c is the child in a trio, p1 is Parent 1 and p2 is Parent 2.)

##Contact
Julie Wertz (julie-wertz@uiowa.edu)
