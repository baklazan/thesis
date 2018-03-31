from alphabet import *
import numpy as np
import sys

class interesting_base:
  def __init__(self, id, reference, correct_value = None):
    self.id = id
    self.reference_value = reference[id]
    if correct_value == None:
      self.real_value = self.reference_value
    else:
      self.real_value = correct_value
    self.log_probability = np.zeros(len(alphabet))
    self.normalized_probability = None
  
  def get_normalized_probability(self):
    if self.normalized_probability is None:
      shifted = (self.log_probability - np.max(self.log_probability))/5
      p = np.exp(shifted)
      sum = np.sum(p)
      self.normalized_probability = p / sum
    return self.normalized_probability
  
  def print(self):
    bullet = ' '
    info = '   '
    if self.real_value != self.reference_value:
      bullet = '*'
      info = "<-{}".format(self.real_value)
    sys.stdout.write("{} {}:\t{}{}".format(bullet, self.id, self.reference_value, info))
    probs = self.get_normalized_probability()
    for i, p in enumerate(probs):
      decoration = " "
      if alphabet[i] == self.reference_value:
        decoration = "*"
      elif alphabet[i] == self.real_value:
        decoration = "-"
      sys.stdout.write("\t{:8.4f}[{:9.4f}]{}".format(p, self.log_probability[i], decoration))
    sys.stdout.write("\n")

  def output(self, f):
    f.write("{} {} {}\n".format(self.real_value, self.reference_value, "\t".join(map(str, self.get_normalized_probability()))))

  def __lt__(self, other):
    return self.id < other.id

def load_interesting_bases(filename, reference):
  result = []
  with open(filename, "r") as f:
    for l in f:
      tokens = l.split()
      result.append(interesting_base(int(tokens[0]), reference, tokens[1]))
  return result
