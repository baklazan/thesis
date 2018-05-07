import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse
from alphabet import *
import copy

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", help="output file for plot", default="plot.png")
parser.add_argument("-s", "--scale", help="portion of plot to be shown", default="0.2")
parser.add_argument("input", help="input file(s)", nargs="+")
args = parser.parse_args()

def plot_rocs(scores, labels):
  plt.xlabel("FAR")
  plt.ylabel("FRR")
  plt.xlim((0, float(args.scale)))
  plt.ylim((0, float(args.scale)))
  points = [(score, labels[i]) for i, score in enumerate(scores)]
  points.sort(reverse=True)
  total_positive = 0
  total_negative = 0
  for label in labels:
    if label == 'positive':
      total_positive += 1
    else:
      total_negative += 1
  x = [0.0]
  y = [1.0]
  best_frr = 1
  false_accepted = 0
  false_rejected = total_positive
  for i, p in enumerate(points):
    if p[1] == 'positive':
      false_rejected -= 1
    else:
      false_accepted += 1
    if i+1 < len(points) and points[i+1][0] == points[i][0]:
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


plt.suptitle("ROC")

for filename in args.input:
  scores = []
  labels = []
  with open(filename, "r") as f:
    for l in f:
      tokens = l.split()
      real_value = tokens[0]
      reference_value = tokens[1]
      p = list(map(float,tokens[2:]))
      is_good = False
      for pr in p:
        is_good = is_good or pr != 0.25
      if is_good:
        score = 1 - p[inv_alphabet[reference_value]]
        label = 'positive' if real_value != reference_value else 'negative'
        scores.append(score)
        labels.append(label)
  plot_rocs(scores, labels)

plt.savefig(args.output, dpi=300)