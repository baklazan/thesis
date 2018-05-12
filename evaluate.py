import matplotlib.pyplot as plt
import numpy as np
import argparse
from alphabet import *
from genome import *
from itertools import chain
import os
import sys
import random

parser = argparse.ArgumentParser()
parser.add_argument("plot_type", help="type of plot to produce")
parser.add_argument("dataset", help="dataset directory")
parser.add_argument("results", help="directories with experiment results", nargs="+")
parser.add_argument("-l", "--labels", help="labels for result files, must come after results", nargs="+")
parser.add_argument("-o", "--output", help="also save plot to file")
args = parser.parse_args()

if args.labels and len(args.labels) != len(args.results):
  print('number of labels must equal number of results')
  sys.exit(1)


def intersection(lists):
  result = []
  for l in lists:
    l.sort()
  indices = [0 for l in lists]
  done = False
  while not done:
    mini = None
    mini_id = None
    good = True
    for i, l in enumerate(lists):
      if indices[i] >= len(l):
        done = True
        break
      if mini != None and l[indices[i]] != mini:
        good = False
      if mini == None or l[indices[i]] < mini:
        mini = l[indices[i]]
        mini_id = i
    if done:
      break
    if good:
      result.append(mini)
    indices[mini_id] += 1
  return result

def plot_rocs(scores, label):
  scores = list(chain.from_iterable(scores))
  scores.sort(reverse=True)
  total_positive = 0
  total_negative = 0
  for score, l in scores:
    if l == 0:
      total_positive += 1
    else:
      total_negative += 1
  x = [0.0]
  y = [0.0]
  best_tpr = 0
  true_positive = 0
  false_positive = 0
  for i, p in enumerate(scores):
    if p[1] == 0:
      true_positive += 1
    else:
      false_positive += 1
    if i+1 < len(scores) and scores[i+1][0] == scores[i][0]:
      continue
    tpr = true_positive / total_positive
    fpr = false_positive / total_negative
    if tpr > best_tpr:
      x.append(fpr)
      y.append(tpr)
      best_tpr = tpr
  x.append(1.0)
  y.append(1.0)
  if len(x) > 100:
    plt.plot(x, y, label=label)
  else:
    plt.plot(x, y, linewidth=0.2, markersize=2, marker='o', label=label)

def decorate_rocs():
  plt.suptitle("ROC")
  plt.xlabel("False positive rate")
  plt.ylabel("True positive rate")
  plt.xlim((0, 1))
  plt.ylim((0, 1))
  plt.legend(loc='lower right')

def plot_distribution(values, label):
  hist, edges = np.histogram(values, bins='auto')
  hist = hist / len(values)
  y = [0, 0]
  x = [0, edges[0]]
  for i, val in enumerate(hist):
    hist[i] = val / (edges[i+1] - edges[i])
    y.append(hist[i])
    x.append(edges[i])
    y.append(hist[i])
    x.append(edges[i+1])
  y.append(0)
  x.append(edges[-1])
  y.append(0)
  x.append(1)
  plt.plot(x, y, label=label)


def plot_score_dist(scores, label):
  scores = list(chain.from_iterable(scores))
  positive = [score for score, type in scores if type == 0]
  negative = [score for score, type in scores if type > 0]
  plot_distribution(positive, "SNP")
  plot_distribution(negative, "nonSNP")

def plot_proximity_to_score_dist(scores, label, classes):
  plt.xlim((0, 1))
  scores = list(chain.from_iterable(scores))
  for c in classes:
    sc = [score for score, dist in filter(lambda x: x[1] >= c[0] and x[1] < c[1], scores)]
    plot_distribution(sc, c[2])

def plot_top_X_success(all_scores, label):
  ratios = []
  for scores in all_scores:
    scores.sort(reverse=True, key=lambda x: (x[0], random.random()))
    SNP_count = len([None for x in scores if x[1] == 0])
    found_SNP_count = len([None for x in scores[:SNP_count] if x[1] == 0])
    if SNP_count > 0:
      ratios.append(found_SNP_count / SNP_count)
  plot_distribution(ratios, label)



def decorate_top_X_success():
  plt.suptitle("Identification success distribution")
  plt.xlim((0, 1))
  plt.ylim(ymin=0)
  plt.xlabel("score")
  plt.ylabel("normalized frequency")
  plt.legend(loc='upper right')


def decorate_score_dist():
  plt.suptitle("Score distribution")
  plt.xlim((0, 1))
  plt.ylim(ymin=0)
  plt.xlabel("score")
  plt.ylabel("normalized frequency")
  plt.legend(loc='upper center')


outfiles = []
for experiment in args.results:
  outfiles.append(list(filter(lambda x: x[-4:] == '.txt', os.listdir(experiment))))

outfiles = intersection(outfiles)

experiments = []
for eid, experiment in enumerate(args.results):
  all_scores = []
  for outfile in outfiles:
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
    nearest_SNP = [float('inf') for c in reference]
    count = float('inf')
    for i in range(len(nearest_SNP)):
      if is_SNP[i]:
        count = 0
      nearest_SNP[i] = min(count, nearest_SNP[i])
      count += 1
    count = float('inf')
    for i in range(len(nearest_SNP)-1, -1, -1):
      if is_SNP[i]:
        count = 0
      nearest_SNP[i] = min(count, nearest_SNP[i])
      count += 1

    scores = []
    with open(os.path.join(experiment, outfile)) as f:
      for l in f:
        tokens = l.split()
        id = int(tokens[0])
        SNP_score = 1 - float(tokens[inv_alphabet[reference[id]]+1])
        scores.append((SNP_score, nearest_SNP[id]))
    all_scores.append(scores)

  label = experiment
  if args.labels:
    label = args.labels[eid]
  experiments.append((all_scores, label))

plot_function = None
decoration_function = None

if args.plot_type == 'ROC':
  plot_function = plot_rocs
  decoration_function = decorate_rocs
elif args.plot_type == 'score_dist':
  plot_function = plot_score_dist
  decoration_function = decorate_score_dist
elif args.plot_type == 'proximity':
  plot_function = lambda x, y: plot_proximity_to_score_dist(x, y, ((0, 1, "SNP"), (1, 2, "next to SNP"), (10, float('inf'), "≥10 bases from SNP"),))
  decoration_function = decorate_score_dist
elif args.plot_type == 'top':
  plot_function = plot_top_X_success
  decoration_function = decorate_top_X_success
else:
  print('unknown metric: {}'.format(args.plot_type))
  sys.exit(1)

for scores, label in experiments:
  plot_function(scores, label)
decoration_function()

if args.output:
  plt.savefig(args.output, dpi=300)
plt.show()
