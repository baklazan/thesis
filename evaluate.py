import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
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
parser.add_argument("-x", "--xlabel", help="label for x axis")
parser.add_argument("-s", "--suptitle", help="suptitle")
parser.add_argument("-z", "--zoom", help="part of ROC curve to be ploted")
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

line_style_cycle = ['-', '--', '-.', ':']
line_style_idx = 0
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
  global line_style_idx
  if len(x) > 100:
    plt.plot(x, y, label=label, linestyle=line_style_cycle[line_style_idx])
  else:
    plt.plot(x, y, linewidth=0.2, markersize=2, marker='o', label=label, linestyle=line_style_cycle[line_style_idx])
  line_style_idx += 1
  line_style_idx %= len(line_style_cycle)

def decorate_rocs():
  plt.suptitle(args.suptitle if args.suptitle else "ROC")
  plt.xlabel(args.xlabel if args.xlabel else "FP")
  plt.ylabel("TP")
  plt.xlim((0, float(args.zoom) if args.zoom else 1))
  plt.ylim((0, 1))
  plt.legend(loc='lower right')



def plot_distribution(values, label, log=False, bins='auto', normed=False):
  color = next(plt.gca()._get_lines.prop_cycler)['color']
  plt.hist(values, label=label, color=to_rgba(color, 0.5), bins=bins, log=log, normed=normed)

def plot_score_dist(scores, label):
  scores = list(chain.from_iterable(scores))
  positive = [score for score, type in scores if type == 0]
  negative = [score for score, type in scores if type > 0]
  print('SNPs: {}'.format(len(positive)))
  print('nonSNPs: {}'.format(len(negative)))
  bins=int(min(len(positive), len(negative)) ** (1/3))
  plot_distribution(positive, "SNP", log=True, bins=bins, normed=True)
  plot_distribution(negative, "bez SNPu", log=True, bins=bins, normed=True)

def plot_proximity_to_score_dist(scores, label, classes):
  scores = list(chain.from_iterable(scores))
  scs = []
  minlen = float('inf')
  for c in classes:
    sc = [score for score, dist in filter(lambda x: x[1] >= c[0] and x[1] < c[1], scores)]
    scs.append(sc)
    minlen = min(minlen, len(sc))

  bins = int(minlen ** (1/3))
  for i, c in enumerate(classes):
    sc = scs[i]
    print('{}: {}'.format(c, len(sc)))
    plot_distribution(sc, c[2], log=True, bins=bins, normed=True)

def plot_top_X_success(all_scores, label, bins='auto'):
  ratios = []
  found_SNP_sum = 0
  SNP_count_sum = 0
  for scores in all_scores:
    scores.sort(reverse=True, key=lambda x: (x[0], random.random()))
    SNP_count = len([None for x in scores if x[1] == 0])
    found_SNP_count = len([None for x in scores[:SNP_count] if x[1] == 0])
    SNP_count_sum += SNP_count
    found_SNP_sum += found_SNP_count
    if SNP_count > 0:
      ratios.append(found_SNP_count / SNP_count)
  plot_distribution(ratios, '{} ({:0.1f} %)'.format(label, found_SNP_sum/SNP_count_sum*100), bins=bins)
  print('{}: {} %'.format(label, found_SNP_sum/SNP_count_sum*100))



def decorate_top_X_success():
  plt.suptitle(args.suptitle if args.suptitle else "Distribúcia úspešnosti identifikácie")
  plt.xlim((0, 1))
  plt.ylim(ymin=0)
  plt.xlabel(args.xlabel if args.xlabel else "skóre")
  plt.ylabel("početnosť")
  plt.legend(loc='lower right')


def decorate_score_dist():
  plt.suptitle(args.suptitle if args.suptitle else "Distribúcia skóre")
  plt.xlim((0, 1))
  plt.xlabel(args.xlabel if args.xlabel else "skóre")
  plt.ylabel("normalizovaná početnosť")
  plt.legend(loc='upper center')


outfiles = []
for experiment in args.results:
  outfiles.append(list(filter(lambda x: x[-4:] == '.txt' and os.stat(os.path.join(experiment, x)).st_size > 0, os.listdir(experiment))))
  emptyfiles = list(filter(lambda x: x[-4:] == '.txt' and os.stat(os.path.join(experiment, x)).st_size == 0, os.listdir(experiment)))
  print('{}: {} files are empty'.format(experiment, len(emptyfiles)))


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
  plot_function = lambda x, y: plot_proximity_to_score_dist(x, y, ((0, 1, "SNP"), (1, 2, "vedľa SNPu"), (10, float('inf'), "≥10 báz od SNPu"),))
  decoration_function = decorate_score_dist
elif args.plot_type == 'top':
  print(min([len(e[0]) for e in experiments]))
  plot_function = lambda x, y: plot_top_X_success(x, y, bins=[i*0.04 for i in range(26)])
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
