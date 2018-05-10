import matplotlib.pyplot as plt
import numpy as np
import argparse
from alphabet import *
from genome import *
from itertools import chain
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", help="also save plot to file")
parser.add_argument("plot_type", help="type of plot to produce")
parser.add_argument("dataset", help="dataset directory")
parser.add_argument("results", help="directory with experiment results")
args = parser.parse_args()


def plot_rocs(scores):
  plt.suptitle("ROC")
  plt.xlabel("FAR")
  plt.ylabel("FRR")
  plt.xlim((0, 1))
  plt.ylim((0, 1))
  scores.sort(reverse=True)
  total_positive = 0
  total_negative = 0
  for score, label in scores:
    if label:
      total_positive += 1
    else:
      total_negative += 1
  x = [0.0]
  y = [1.0]
  best_frr = 1
  false_accepted = 0
  false_rejected = total_positive
  for i, p in enumerate(scores):
    if p[1]:
      false_rejected -= 1
    else:
      false_accepted += 1
    if i+1 < len(scores) and scores[i+1][0] == scores[i][0]:
      continue
    frr = false_rejected / total_positive
    far = false_accepted / total_negative
    if frr < best_frr:
      x.append(far)
      y.append(frr)
      best_frr = frr
  x.append(1.0)
  y.append(0.0)
  if len(x) > 100:
    plt.plot(x, y)
  else:
    plt.plot(x, y, linewidth=0.2, markersize=2, marker='o')

def plot_distribution(values):
  hist, edges = np.histogram(values, bins=15)
  hist = hist / len(values)
  x = [(edges[i] + edges[i+1])/2 for i in range(len(edges)-1)]
  plt.plot(x, hist)

def plot_score_dist(scores):
  positive = [score for score, type in scores if type]
  negative = [score for score, type in scores if not type]
  plot_distribution(positive)
  plot_distribution(negative)

all_scores = []
for outfile in os.listdir(args.results):
  if outfile[-4:] != '.txt':
    continue
  reffile = os.path.join(args.dataset, 'reads', outfile[:-4], 'reference.fasta')
  reference = load_fasta(reffile)[0].bases
  snpfile = os.path.join(args.dataset, 'reads', outfile[:-4], 'changes.txt')
  is_SNP = [False for c in reference]
  true_value = reference[:]
  with open(snpfile, 'r') as f:
    for l in f:
      id, val = l.split()
      id = int(id)
      true_value[id] = val
      is_SNP[id] = True

  scores = []
  with open(os.path.join(args.results, outfile)) as f:
    for l in f:
      tokens = l.split()
      id = int(tokens[0])
      SNP_score = 1 - float(tokens[inv_alphabet[reference[id]]+1])
      scores.append((SNP_score, is_SNP[id]))
  all_scores.append(scores)

if args.plot_type == 'ROC':
  scores = list(chain.from_iterable(all_scores))
  plot_rocs(scores)
elif args.plot_type == 'score_dist':
  scores = list(chain.from_iterable(all_scores))
  plot_score_dist(scores)
else:
  print('unknown metric: {}'.format(args.plot_type))
  sys.exit(1)


if args.output:
  plt.savefig(args.output, dpi=300)
plt.show()
