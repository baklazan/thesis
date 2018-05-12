import h5py
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
import argparse
from genome import *
from read import *
from bc_read import *
from alphabet import *
from kmer_model import *
from base import *
from model import *
import yaml
import os

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", help="output file for a posteriori probabilities (default: stdout)")
parser.add_argument("-i", "--independent", help="treat each read independently", action='store_true')
parser.add_argument("-k", "--kmer_model", help="file with kmer model to use", default="tombo")
parser.add_argument("-g", "--group_name", help="name of group in fast5 files containing basecall info", default="Analyses/Basecall_1D_000")
parser.add_argument("-c", "--configuration", help="config file with reresquiggle parameters", default="default.yaml")
parser.add_argument("reference", help="reference fasta file")
parser.add_argument("read_basedir", help="base directory of resquiggled fast5 files")

args = parser.parse_args()


kmer_model = None
if args.kmer_model == "tombo":
  kmer_model = tombo_kmer_model("data/tombo.DNA.model")
elif args.kmer_model == "picoamp":
  kmer_model = picoamp_kmer_model("data/6mer_model.txt")
elif args.kmer_model == "klebs":
  kmer_model = tombo_kmer_model("data/klebs2of6.DNA.model")
else:
  raise ValueError("Unknown kmer model: {}".format(args.kmer_model))

with open(args.configuration, "r") as f:
  config = yaml.load(f)

load_read = None
if config['model'] == 'basecall':
  model = basecall_model(args.reference, config = config)
  load_read = lambda read_file, kmer_model, group_name : basecalled_read(filename=read_file, kmer_model = kmer_model, basecall_group=group_name)
elif config['model'] == 'windowed':
  model = window_model(kmer_model, config = config)
  load_read = lambda read_file, kmer_model, group_name : resquiggled_read(filename=read_file, kmer_model = kmer_model)
elif config['model'] == 'moving':
  model = moving_window_model(kmer_model, config = config)
  load_read = lambda read_file, kmer_model, group_name: resquiggled_read(filename=read_file, kmer_model=kmer_model)
else:
  print('Unknown model: {}'.format(config['model']))
  sys.exit(1)

try:
  reference = load_fasta(args.reference)[0].bases
except FileNotFoundError:
  print("failed to process: reference {} doesn't exist".format(args.reference))
  exit(1)
interesting = [interesting_base(i, reference) for i, _ in enumerate(reference)]


read_files = [os.path.join(args.read_basedir, file) for file in os.listdir(args.read_basedir) if not os.path.isdir(os.path.join(args.read_basedir, file))]
read_files = filter(lambda x : x[-6:] == ".fast5", read_files)

outfile = sys.stdout
if args.output:
  outfile = open(args.output, 'w')

reads = []
for read_file in read_files:
  try:
    read = load_read(read_file, kmer_model, args.group_name)
  except KeyError:
    print('failed to process read {}'.format(read_file))
    continue
  read.fix_start_in_reference(reference)
  if config['tweak_normalization']:
    if read.strand == '-':
      read.tweak_normalization(reverse_complement(reference), kmer_model)
    else:
      read.tweak_normalization(reference, kmer_model)
  print("[{}, {}){}".format(read.start_in_reference, read.end_in_reference, read.strand))

  model.update_probabilities(reference, read, interesting)
  reads.append(read)
  if args.independent:
    for base in interesting:
      base.output(outfile)
      base.clear_probabilities()

if not args.independent:
  for base in interesting:
    base.output(outfile)

if args.output:
  outfile.close()