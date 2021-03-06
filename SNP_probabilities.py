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
parser.add_argument("-o", "--output", help="output file for a posteriori probabilities")
parser.add_argument("-l", "--log_output", help="file for detailed output (including conditional probabilities and base indices)")
parser.add_argument("-p", "--print_worst", help="print worst n candidates")
parser.add_argument("-s", "--show_worst", help="number of worst positions to show (only works with --print_best)")
parser.add_argument("-i", "--independent", help="treat each read independently", action='store_true')
parser.add_argument("-k", "--kmer_model", help="file with kmer model to use", default="tombo")
parser.add_argument("-g", "--group_name", help="name of group in fast5 files containing basecall info", default="Analyses/Basecall_1D_000")
parser.add_argument("-b", "--basecall_only", help="only use basecalled sequence to determine SNPs (ignore signal)", action='store_true')
parser.add_argument("-c", "--configuration", help="config file with reresquiggle parameters", default="default.yaml")
parser.add_argument("-f", "--full_interval", help="compute scores for whole reference instead of just interesting bases", action="store_true")
parser.add_argument("-a", "--around_changes", help="compute scores for all bases between first and last changed base", action="store_true")
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
elif args.kmer_model == "klebs":
  kmer_model = tombo_kmer_model("data/klebs2of6.DNA.model")
else:
  raise ValueError("Unknown kmer model: {}".format(args.kmer_model))

with open(args.configuration, "r") as f:
  config = yaml.load(f)

load_read = None
if args.basecall_only:
  model = basecall_model(args.reference, config)
  load_read = lambda read_file, kmer_model, group_name : basecalled_read(filename=read_file, kmer_model = kmer_model, basecall_group=group_name)

else:
  model = window_model(kmer_model, config = config)
  load_read = lambda read_file, kmer_model, group_name : resquiggled_read(filename=read_file, kmer_model = kmer_model)

reference = load_fasta(args.reference)[0].bases

interesting = load_interesting_bases(args.interesting, reference)
if args.full_interval and args.around_changes:
  print("Can't use both --aroung_changes and --full_interval")
  exit(1)
if args.full_interval:
  new_interesting = [None for b in reference]
  for b in interesting:
    new_interesting[b.id] = b
  interesting = new_interesting
  for i in range(len(reference)):
    if interesting[i] == None:
      interesting[i] = interesting_base(i, reference)
elif args.around_changes:
  changed = list(filter(lambda b: b.real_value != b.reference_value, interesting))
  begin = max(0, min(changed).id - 20)
  end = min(len(reference), max(changed).id + 1 + 20)
  interesting = [None for i in range(begin, end)]
  for b in changed:
    interesting[b.id - begin] = b
  for i in range(begin, end):
    if interesting[i - begin] == None:
      interesting[i - begin] = interesting_base(i, reference)


read_files = [os.path.join(args.read_basedir, file) for file in os.listdir(args.read_basedir) if not os.path.isdir(os.path.join(args.read_basedir, file))]
read_files = filter(lambda x : x[-6:] == ".fast5", read_files)

if args.output:
  open(args.output, "w").close()
if args.log_output:
  open(args.log_output, "w").close()

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

  model.update_probabilities_full(reference, read, interesting)
  reads.append(read)
  if args.independent:
    if args.output:
      with open(args.output, "a") as f:
        for base in interesting:
          base.output(f)
    if args.log_output:
      with open(args.log_output, "a") as f:
        for base in interesting:
          base.log_output(f)
    for base in interesting:
      base.clear_probabilities()

if not args.independent:
  for base in interesting:
    if base.real_value != base.reference_value:
      base.print()
  print()

  if args.full_interval or args.around_changes:
    x, prob, conf, maxconf = [], [], [], []
    for b in interesting:
      probs = b.get_normalized_probability()
      prob.append(1 - probs[inv_alphabet[b.reference_value]])
      conf.append(b.log_probability[inv_alphabet[b.reference_value]])
      maxconf.append(max(b.log_probability))
      x.append(b.id)


    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)

    plt.plot(x, prob)
    #pl2 = plt.twinx()
    #pl2.plot(x, conf)
    #pl2.plot(x, maxconf, color = 'red')
    for b in interesting:
      if b.real_value != b.reference_value:
        plt.axvspan(b.id-0.5, b.id+0.5, color='green', alpha=0.5)

    from matplotlib.widgets import Slider


    axcolor = 'lightgoldenrodyellow'
    axpos = plt.axes([0.2, 0.1, 0.65, 0.03], facecolor=axcolor)

    spos = Slider(axpos, 'Pos', x[0], x[-1] - 100 + 1)


    def update(val):
      pos = val
      ax.axis([pos, pos + 100, 0, 1])
      fig.canvas.draw_idle()

    spos.on_changed(update)
    plt.show()

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
  if args.log_output:
    with open(args.log_output, "w") as f:
      for base in interesting:
        base.log_output(f)
