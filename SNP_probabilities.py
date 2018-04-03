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
import os

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", help="output file for a posteriori probabilities")
parser.add_argument("-p", "--print_worst", help="print worst n candidates")
parser.add_argument("-s", "--show_worst", help="number of worst positions to show (only works with --print_best)")
parser.add_argument("-i", "--independent", help="treat each read independently", action='store_true')
parser.add_argument("-k", "--kmer_model", help="file with kmer model to use", default="tombo")
parser.add_argument("reference", help="reference fasta file")
parser.add_argument("interesting", help="list of interesting positions in reference")
parser.add_argument("read_basedir", help="base directory of resquiggled fast5 files")

args = parser.parse_args()

if args.independent and args.print_worst:
  print("--print_worst doesn't work with --independent")
  exit(1)

if args.show_worst and not args.print_worst:
  print("--show_worst only works with --print_worst")
  exit(1)

kmer_model = None
if args.kmer_model == "tombo":
  kmer_model = tombo_kmer_model("data/tombo.DNA.model")
elif args.kmer_model == "picoamp":
  kmer_model = picoamp_kmer_model("data/6mer_model.txt")
else:
  raise ValueError("Unknown kmer model: {}".format(args.kmer_model))
model = window_model(kmer_model, op = max_operation(), buffer_size=8, min_event_length=2, window_size=23, penalty=0.5)
reference = load_fasta(args.reference)[0].bases
interesting = load_interesting_bases(args.interesting, reference)

read_files = [os.path.join(args.read_basedir, file) for file in os.listdir(args.read_basedir) if not os.path.isdir(os.path.join(args.read_basedir, file))]

open(args.output, "w").close()

reads = []
for read_file in read_files:
  try:
    read = resquiggled_read(read_file, kmer_model)
  except KeyError:
    continue
  read.fix_start_in_reference(reference)
  if read.strand == '-':
    read.tweak_normalization(reverse_complement(reference), kmer_model)
  else:
    read.tweak_normalization(reference, kmer_model)
  print("[{}, {})".format(read.start_in_reference, read.end_in_reference))

  model.update_probabilities_c(reference, read, interesting)
  reads.append(read)
  if args.independent:
    if args.output:
      with open(args.output, "a") as f:
        for base in interesting:
          base.output(f)
    for base in interesting:
      base.clear_probabilities()

if not args.independent:
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

  if args.output:
    with open(args.output, "w") as f:
      for base in interesting:
        base.output(f)
