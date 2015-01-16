# Exome-VCF-processing
vcf_compare.py, vcf_compare_inheritance.py

1/16/15

##Purpose
Compares variants between any number of VCF files. Variant comparison files,
QC plots, and statistics are produced.

vcf_compare_inheritance.py: Additionaly analyzes inheritance.
(Used with VCF files obtained from a parent-offspring trio.)

##Requirements
Python 2.x (https://www.python.org/)

PyVCF (https://github.com/jamescasbon/PyVCF)

matplotlib (http://matplotlib.org/)

##Usage
###vcf_compare.py:

####Arguments: 
vcf_file_names: Names of VCF files to be processed (required)

####Options:
-h, --help: Show help message and exit.

--create_qc_plots: Whether to output QC plots. Takes no arguments.
(Default: Does not output plots.)

--out_stats: Name of HTML table output file containing overall statistics.
(Default: vcf_stats.html)

--depth_plot: Name of output file containing read depth QC plot.
(--create_qc_plots must be specified.) (Default: dp_plot.png)

--qd_plot: Name of output file containing quality by depth QC plot.
(--create_qc_plots must be specified.) (Default: qd_plot.png)

--qual_plot: Name of output file containing quality QC plot.
(--create_qc_plots must be specified.) (Default: qual_plot.png)

--out_all_variants: Name of output file containing all high-quality variants.
(Default: all_variants.txt)

--out_absent: Name of output file containing variants such that at least one
file has high-quality variant and another file does not have variant at CHROM
and POS (Default: conflicts_absent.txt)

--out_multiple: Name of output file containing variants such that at least one
file has high-quality variant and another file has different REF/ALT/GT at
CHROM and POS. Excludes indels. (Default: conflicts_multiple.txt)

--min_qd: Minimum quality by depth for inclusion in variant output files.
(Default: 3)

--min_qual: Minimum quality score for inclusion in variant output files.
(Default: 80)

--min_dp: Minimum read depth for inclusion in variant output files.
(Default: 15)

####Example commands:

python vcf_compare.py --create_qc_plots a.vcf b.vcf c.vcf --qd_plot d.png
--depth_plot e.png --qual_plot f.png

python vcf_compare.py a.vcf b.vcf --out_stats c.html --out_all_variants
d.txt --out_absent e.txt --out_multiple f.txt --min_qd 5 --min_qual 100
--min_dp 20
 
###vcf_compare_inheritance.py:

####Arguments: 

-p1, --parent_1_file_names: VCF filenames corresponding to Parent 1
(required)

-p2, --parent_2_file_names: VCF filenames corresponding to Parent 2
(required)

-c, --child_file_names: VCF filenames corresponding to child (required)

####Options:
Same as above, excluding out_absent and out_multiple.

####Example commands: 
python vcf_compare_inheritance.py --create_qc_plots -p1 a.vcf b.vcf
-p2 c.vcf d.vcf -c e.vcf f.vcf --qd_plot g.png --depth_plot h.png
--qual_plot i.png

python vcf_compare_inheritance.py -p1 a.vcf -p2 b.vcf -c c.vcf
--out_stats d.html --out_all_variants e.txt --min_qd 5 --min_qual 100
--min_dp 20

##Contact
Julie Wertz (julie-wertz@uiowa.edu)
