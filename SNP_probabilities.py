import h5py
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
import argparse
from genome import *
from read import *
from alphabet import *
from kmer_model import *
from base import *
from model import *

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", help="output file for a posteriori probabilities")
parser.add_argument("-p", "--print_worst", help="print worst n candidates")
parser.add_argument("-k", "--kmer_model", help="file with kmer model to use", default="tombo")
parser.add_argument("-s", "--show_worst", help="number of worst positions to show (only works with --print_best)", default="tombo")
parser.add_argument("reference", help="reference fasta file")
parser.add_argument("interesting", help="list of interesting positions in reference")
parser.add_argument("read", help="resquiggled fast5 file(s)", nargs="+")
args = parser.parse_args()


kmer_model = None
if args.kmer_model == "tombo":
  kmer_model = tombo_kmer_model("data/tombo.DNA.model")
elif args.kmer_model == "picoamp":
  kmer_model = picoamp_kmer_model("data/6mer_model.txt")
else:
  raise ValueError("Unknown kmer model: {}".format(args.kmer_model))
model = window_model(kmer_model, op = max_operation(), buffer_size=4, min_event_length=2, max_event_length=60, window_size=13)
reference = load_fasta(args.reference)[0].bases
interesting = load_interesting_bases(args.interesting, reference)

reads = []
for read_file in args.read:
  read = resquiggled_read(read_file, kmer_model)
  print("[{}, {})".format(read.start_in_reference, read.end_in_reference))
  read.tweak_normalization(reference, kmer_model)
  model.update_probabilities(reference, read, interesting)
  reads.append(read)

for base in interesting:
  if base.real_value != base.reference_value:
    base.print()
print()

if args.print_worst:
  by_SNP_prob = []
  for base in interesting:
    probs = base.get_normalized_probability()
    by_SNP_prob.append((1 - probs[inv_alphabet[base.real_value]], base))
  by_SNP_prob.sort(reverse=True)
  for p, base in by_SNP_prob[:int(args.print_worst)]:
    base.print()
  if args.show_worst:
    for p, base in by_SNP_prob[:int(args.show_worst)]:
      for read in reads:
        model.show_base(reference, read, base)

count_good = 0
count_bad = 0
for base in interesting:
  p = base.log_probability
  if p[inv_alphabet[base.reference_value]] == np.max(p):
    count_good += 1
  else:
    count_bad += 1

print("good: {} / {} [{} %]".format(count_good, count_good+count_bad, count_good/(count_good + count_bad)*100))

if args.output:
  with open(args.output, "w") as f:
    for base in interesting:
      base.output(f)