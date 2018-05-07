import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse
from alphabet import *
import copy

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", help="output file for plot", default="plot.png")
parser.add_argument("input", help="input file(s)", nargs="+")
args = parser.parse_args()

epsilon = 1e-10

def confusion_matrix(scores, labels, weights):
  result = np.zeros((len(alphabet), len(alphabet)))
  for i, l in enumerate(labels):
    id = inv_alphabet[l]
    result[id, np.argmax(scores[i] * weights)] += 1
  return result

def all_trees(n):
  if n != 4:
    raise ValueError("number of classes isn't 4")
  return [[(0, 1), (0, 2), (0, 3)],
          [(0, 1), (0, 2), (1, 3)],
          [(0, 1), (0, 2), (2, 3)],
          [(0, 1), (0, 3), (1, 2)],
          [(0, 1), (0, 3), (2, 3)],
          [(0, 1), (1, 2), (1, 3)],
          [(0, 1), (1, 2), (2, 3)],
          [(0, 1), (1, 3), (2, 3)],
          [(0, 2), (0, 3), (1, 2)],
          [(0, 2), (0, 3), (1, 3)],
          [(0, 2), (1, 2), (1, 3)],
          [(0, 2), (1, 2), (2, 3)],
          [(0, 2), (2, 3), (1, 3)],
          [(0, 3), (1, 3), (1, 2)],
          [(0, 3), (2, 3), (1, 2)],
          [(0, 3), (1, 3), (2, 3)]]

def generate_weights_data(tree, tresholds, current_weights, current_edge, result):
  if current_edge == len(tree):
    result.append(copy.deepcopy(current_weights))
    return
  u = tree[current_edge][0]
  v = tree[current_edge][1]
  for t in tresholds[u][v]:
    if current_weights[u] != None:
      current_weights[v] = current_weights[u] * (1-t) / t
      generate_weights_data(tree, tresholds, current_weights, current_edge+1, result)
      current_weights[v] = None
    elif current_weights[v] != None:
      current_weights[u] = current_weights[v] * t / (1-t)
      generate_weights_data(tree, tresholds, current_weights, current_edge+1, result)
      current_weights[u] = None
    else:
      current_weights[u] = t
      current_weights[v] = 1-t
      generate_weights_data(tree, tresholds, current_weights, current_edge+1, result)
      current_weights[u] = None
      current_weights[v] = None

resolution = 40
def generate_weights_uniform(current_weights, index, result, sum = 0):
  if index+1 == len(current_weights):
    current_weights[index] = (resolution - sum) / resolution
    result.append(copy.deepcopy(current_weights))
    return
  for x in range(1, resolution):
    if x + sum >= resolution:
      break
    w = x / resolution
    current_weights[index] = w
    generate_weights_uniform(current_weights, index+1, result, sum + x)

def roc_points(scores, labels):
  result = []
  if len(labels) < 5:
    tresholds = [[[] for i in alphabet] for j in alphabet]
    for i in range(len(alphabet)):
      for j in range(i+1, len(alphabet)):
        for k in scores:
          if k[i] > 0 and k[j] > 0:
            t = k[j] / (k[i] + k[j])
            tresholds[i][j].append(t - epsilon)
            tresholds[i][j].append(t + epsilon)
    trees = all_trees(len(alphabet))
    for tree in trees:
      weights = []
      generate_weights_data(tree, tresholds, [None for i in range(len(alphabet))], 0, weights)
      for w in weights:
        result.append(confusion_matrix(scores, labels, w))
  else:
    weights = []
    generate_weights_uniform([None for i in range(len(alphabet))], 0, weights)
    for w in weights:
      result.append(confusion_matrix(scores, labels, w))
  return result

def plot_rocs(scores, labels, plots):
  points = roc_points(scores, labels)
  for i, plot in enumerate(plots):
    plot.set_xlabel("FAR")
    plot.set_ylabel("FRR")
    plot.set_title(alphabet[i])
    plot.set_xlim((0, 0.2))
    plot.set_ylim((0, 0.2))
    my_points = []
    for p in points:
      p_was_me = np.sum(p[i])
      frr = (p_was_me - p[i][i]) / p_was_me
      p_missclassified_me = np.sum(p[:, i]) - p[i][i]
      p_wasnt_me = np.sum(p) - p_was_me
      far = p_missclassified_me / p_wasnt_me
      my_points.append((far, frr))
    my_points.sort()
    x = [0.0]
    y = [1.0]
    best_frr = 1
    for p in my_points:
      if p[1] < best_frr:
        x.append(p[0])
        y.append(p[1])
        best_frr = p[1]
    x.append(1.0)
    y.append(0.0)
    plot.plot(x, y)

rows = np.floor(np.sqrt(len(alphabet)))
columns = np.ceil(len(alphabet) / rows)
plt.figure()
plt.suptitle("ROCs")
plots = [plt.subplot(rows, columns, i+1) for i in range(len(alphabet))]
for filename in args.input:
  scores = []
  labels = []
  with open(filename, "r") as f:
    for l in f:
      tokens = l.split()
      labels.append(tokens[0])
      scores.append(np.array(list(map(float,tokens[2:]))))
  plot_rocs(scores, labels, plots)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(args.output, dpi=300)