import argparse
from genome import *

parser = argparse.ArgumentParser()
parser.add_argument("reference", help="fasta file with reference genome")
parser.add_argument("interesting", help="list of interesting positions in reference")
parser.add_argument("log", help="file containing log of basecall-based model")


args = parser.parse_args()
reference = load_fasta(args.reference)[0].bases

real_letter = {}
with open(args.interesting, 'r') as f:
  for l in f:
    tokens = l.split()
    real_letter[int(tokens[0])] = tokens[1]

non_SNP_count = 0
SNP_count = 0

SNP_insert = 0
SNP_delete = 0
SNP_false_good = 0
SNP_true_substitute = 0
SNP_false_substitute = 0

non_SNP_insert = 0
non_SNP_delete = 0
non_SNP_good = 0
non_SNP_substitute = 0

with open(args.log, 'r') as f:
  for l in f:
    tokens = l.split()
    id = int(tokens[0])
    is_SNP = real_letter[id] != reference[id]
    if is_SNP:
      SNP_count += 1
    else:
      non_SNP_count += 1
    if tokens[1] == 'I':
      if is_SNP:
        SNP_insert += 1
      else:
        non_SNP_insert += 1
    elif tokens[1] == 'D':
      if is_SNP:
        SNP_delete += 1
      else:
        non_SNP_delete += 1
    elif tokens[1] == 'G':
      if is_SNP:
        SNP_false_good += 1
      else:
        non_SNP_good += 1
    elif tokens[1] == 'S':
      if is_SNP:
        if tokens[2] == real_letter[id]:
          SNP_true_substitute += 1
        else:
          SNP_false_substitute += 1
      else:
        non_SNP_substitute += 1

print('SNP_insert: {}'.format(SNP_insert / SNP_count))
print('SNP_delete: {}'.format(SNP_delete / SNP_count))
print('SNP_false_good: {}'.format(SNP_false_good / SNP_count))
print('SNP_true_substitute: {}'.format(SNP_true_substitute / SNP_count))
print('SNP_false_substitute: {}'.format(SNP_false_substitute / SNP_count))

print('non_SNP_insert: {}'.format(non_SNP_insert / SNP_count))
print('non_SNP_delete: {}'.format(non_SNP_delete / SNP_count))
print('non_SNP_good: {}'.format(non_SNP_good / SNP_count))
print('non_SNP_substitute: {}'.format(non_SNP_substitute / SNP_count))
