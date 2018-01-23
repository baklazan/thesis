import h5py
import numpy as np
import sys
import random
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--changed", help="print list of changed positions into file")
parser.add_argument("-n", "--number", help="number of changes to make")
parser.add_argument("-p", "--percentage", help="percentage of changes to make (doesn't work if --number is used)")
parser.add_argument("reference", help="reference fasta file")
parser.add_argument("output", help="output fasta file")
parser.add_argument("-r", "--reads", help="only change bases in range of these reads", nargs="+")
args = parser.parse_args()


alphabet = ['A', 'C', 'G', 'T']
def something_else(base):
  res = alphabet[random.randrange(len(alphabet))]
  while res == base:
    res = alphabet[random.randrange(len(alphabet))]
  return res

reference = ""
reference_head = []
reference_tail = []
reference_lines = []
with open(args.reference, "r") as f:
  lines = []
  first_started = False
  in_tail = False
  for l in f:
    if l[0] == '>':
      if first_started:
        in_tail = True
        reference_tail.append(l)
      else:
        first_started = True
        reference_head.append(l)
    else:
      if in_tail:
        reference_tail.append(l)
      else:
        lines.append(l.rstrip())
  reference = "".join(lines)

valid_positoins = []
if args.reads:
  valid_ranges = []
  for filename in args.reads:
    with h5py.File(filename) as f:
      alignment_meta = f['Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment'].attrs
      valid_ranges.append((alignment_meta['mapped_start'], alignment_meta['mapped_end']))
  dstate = np.zeros(len(reference)+1)
  for r in valid_ranges:
    dstate[r[0]] += 1
    dstate[r[1]] -= 1
  sum = 0
  for i, d in enumerate(dstate):
    sum += d
    if sum > 0:
      valid_positoins.append(i)
else:
  valid_positoins = list(range(len(reference)))

ref_array = list(reference)

changes = []
if args.number:
  for i in range(int(args.number)):
    index = valid_positoins[random.randrange(len(valid_positoins))]
    changes.append((index, ref_array[index]))
    ref_array[index] = something_else(ref_array[index])
elif args.percentage:
  for index in valid_positoins:
    if random.randrange(100) < args.percentage:
      changes.append((index, ref_array[index]))
      ref_array[index] = something_else(ref_array[index])
else:
  print("Neither --number nor --percentage were specified")
  sys.exit(1)

with open(args.output, "w") as f:
  for l in reference_head:
    f.write(l)
  for i in range(0, len(ref_array), 70):
    f.write("{}\n".format("".join(ref_array[i:i+70])))
  for l in reference_tail:
    f.write(l)

if args.changed:
  with open(args.changed, "w") as f:
    for c in changes:
      f.write("{} {}\n".format(c[0], c[1]))