from kmer_model import *
from base import *
import random

class max_operation:
  def __init__(self):
    pass

  def neutral_element(self):
    return -np.inf

  def combine(self, a, b):
    return max(a, b)

  def start_element(self):
    return 0

  def append_combine(self, a, b):
    return a + b


class log_sum_operation:
  def __init__(self):
    pass

  def neutral_element(self):
    return -np.inf

  def combine(self, a, b):
    if np.isinf(a) and np.isinf(b):
      return -np.inf
    if a < b:
      return b + np.log(1 + np.exp(a-b))
    else:
      return a + np.log(1 + np.exp(b-a))

  def start_element(self):
    return 0

  def append_combine(self, a, b):
    return a + b


class window_model:
  def __init__(self, kmer_model = None, window_size = 15, min_event_length = 3, max_event_length = 30, op = None):
    self.window_size = window_size
    self.min_event_length = min_event_length
    self.max_event_length = max_event_length
    self.kmer_model = kmer_model
    self.op = op
    if op == None:
      self.op = max_operation()

  def reresquiggle(self, signal, reference_kmer_ids, start_buffer=0, end_buffer=0):
    print("st: {} en: {}".format(start_buffer, end_buffer))
    dp = np.full((len(signal) + 1, len(reference_kmer_ids) + 1), self.op.neutral_element(), dtype=float)
    for j in range(start_buffer+1):
      dp[j][0] = self.op.start_element()
    for i, kmer_id in enumerate(reference_kmer_ids):
      for j in range(len(signal) + 1):
        accumulator = 0
        for l in range(1, min(self.max_event_length, j) + 1):
          p = self.kmer_model.log_chance_of_signal(kmer_id, signal[j-l])
          accumulator = self.op.append_combine(p, accumulator)
          if l >= self.min_event_length:
            dp[j][i + 1] = self.op.combine(dp[j][i+1], self.op.append_combine(dp[j-l][i], accumulator))
    res = np.max(dp[-end_buffer-1:, -1])
    return res

  def update_base(self, reference, base, signal_values, start_buffer, end_buffer):
    ref_st = base.id - self.window_size // 2 - self.kmer_model.central_position
    ref_en = ref_st + self.window_size + self.kmer_model.k - 1
    if ref_st < 0 or ref_en >= len(reference):
      return
    context = list(reference[ref_st+1:ref_en-1])
    for i, c in enumerate(alphabet):
      context[self.kmer_model.central_position + self.window_size // 2 - 1] = c
      kmer_ids = [kmer_to_id(context[j:j + self.kmer_model.k]) for j in range(len(context) - self.kmer_model.k - 1)]
      base.log_probability[i] += self.reresquiggle(signal_values, kmer_ids, start_buffer, end_buffer)

  def update_probabilities(self, reference, read, interesting_bases):
    for base in interesting_bases:
      id_in_read = base.id - read.start_in_reference
      if id_in_read - self.window_size // 2 < 0 or id_in_read + self.window_size // 2 >= len(read.event_start):
        continue
      start = read.event_start[id_in_read - self.window_size // 2]
      end = read.event_start[id_in_read + self.window_size // 2] + read.event_length[id_in_read + self.window_size // 2]
      signals = read.normalized_signal[start:end]
      start_buffer = read.event_start[id_in_read - self.window_size // 2 + 2] - start
      end_buffer = end - read.event_start[id_in_read + self.window_size // 2 - 1]
      self.update_base(reference, base, signals, start_buffer, end_buffer)

class no_model:
  def __init__(self, kmer_model = None):
    pass

  def update_base(self, base):
    for i in range(len(alphabet)):
      base.log_probability[i] -= random.random()

  def update_probabilities(self, reference, read, interesting_bases):
    for base in interesting_bases:
      self.update_base(base)
