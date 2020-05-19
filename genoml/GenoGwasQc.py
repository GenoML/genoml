import argparse
import sys

#local imports
import genoml.dependencies
from genoml.gwas import qc_steps


parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')
parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
# parser.add_argument('--out', type=str, default='out', help='Prefix for output (including path)')
parser.add_argument('--rare', default=False, action="store_true", help='Pruning toggle for rare variants. If --rare is used, final MAF pruning (0.01) will not be conducted, otherwise, rare variants will be pruned')

args = parser.parse_args()

geno_name = args.geno
rare_flag = args.rare

# INSTANTIATE QC WITH INPUT NAME AND OUTPUT NAME     
qc = qc_steps.QC(geno_name)

# NOW RUN COMMANDS
# FIRST, CLEAR EXISTING LOGFILE
qc.rm_log()

# run het pruning
qc.het_pruning()
qc.call_rate_pruning()
qc.sex_check()
qc.relatedness_pruning()
qc.variant_pruning()
qc.rare_prune()
qc.cleanup()

print("DONE!!!!!!!")
