import h5py
import numpy as np
import sys
import os
import random
import argparse
from alphabet import *
from genome import *

parser = argparse.ArgumentParser()
parser.add_argument('reads_basedir', help='parent of directories containing reads')
parser.add_argument('reference', help='reference fasta file')
parser.add_argument('percentage', help='percentage of changes to make')
args = parser.parse_args()


def something_else(base):
  res = alphabet[random.randrange(len(alphabet))]
  while res == base:
    res = alphabet[random.randrange(len(alphabet))]
  return res

reference = load_fasta(args.reference)[0].bases
read_dirs = [os.path.join(args.reads_basedir, dir) for dir in os.listdir(args.reads_basedir)
             if os.path.isdir(os.path.join(args.reads_basedir, dir))]
for dir in read_dirs:
  read_file = os.path.join(dir, 'read.fast5')
  start, end = None, None
  try:
    with h5py.File(read_file) as f:
      alignment_meta = f['Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment'].attrs
      start, end = alignment_meta['mapped_start'], alignment_meta['mapped_end']
  except KeyError:
    print('failed to process read {}'.format(dir))
    continue
  start = max(0, start - 50)
  end = min(len(reference), end + 50)
  relevant_reference = reference[start:end]
  changes = []
  for i, c in enumerate(relevant_reference):
    if random.random()*100 < float(args.percentage):
      relevant_reference[i] = something_else(c)
      changes.append((i, relevant_reference[i]))
  with open(os.path.join(dir, 'reference.fasta'), 'w') as f:
    f.write('>Interval [{}, {}) from original reference\n'.format(start, end))
    for line_start in range(0, len(relevant_reference), 80):
      line_end = min(line_start+80, len(relevant_reference))
      f.write('{}\n'.format(''.join(relevant_reference[line_start:line_end])))
  with open(os.path.join(dir, 'changes.txt'), 'w') as f:
    for pos, val in changes:
      f.write('{}\t{}\n'.format(pos, val))
