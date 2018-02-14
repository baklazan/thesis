from kmer_model import *
from base import *

class window_model:
  def __init__(self, kmer_model = None, window_size = 11, min_event_length = 2, max_event_length = 30):
    self.window_size = window_size
    self.min_event_length = min_event_length
    self.max_event_length = max_event_length
    self.kmer_model = kmer_model

  def reresquiggle(self, signal, reference_kmer_ids):
    prefix_sums = np.zeros(len(signal) + 1)
    for i, s in enumerate(signal):
      prefix_sums[i + 1] = prefix_sums[i] + signal[i]
    dp = np.full((len(signal) + 1, len(reference_kmer_ids) + 1), -np.inf, dtype=float)
    dp[0][0] = 0
    for i, kmer_id in enumerate(reference_kmer_ids):
      for j in range(len(signal) + 1):
        for l in range(self.min_event_length, min(self.max_event_length, j) + 1):
          average = (prefix_sums[j] - prefix_sums[j - l]) / l
          dp[j][i + 1] = max(dp[j][i + 1], self.kmer_model.log_chance_of_signal(kmer_id, average) + dp[j - l][i])
    return dp[len(signal)][len(reference_kmer_ids)]

  def update_base(self, reference, base, signal_values):
    ref_st = base.id - self.window_size // 2 - self.kmer_model.central_position
    ref_en = ref_st + self.window_size + self.kmer_model.k - 1
    if ref_st < 0 or ref_en >= len(reference):
      return
    context = list(reference[ref_st:ref_en])
    for i, c in enumerate(alphabet):
      context[self.kmer_model.central_position + self.window_size // 2] = c
      kmer_ids = [self.kmer_model.kmer_to_id(context[j:j + self.kmer_model.k]) for j in range(len(context) - self.kmer_model.k - 1)]
      base.log_probability[i] += self.reresquiggle(signal_values, kmer_ids)

  def update_probabilities(self, reference, read, interesting_bases):
    for base in interesting_bases:
      id_in_read = base.id - read.start_in_reference
      if id_in_read - self.window_size // 2 < 0 or id_in_read + self.window_size // 2 >= len(read.event_start):
        continue
      start = read.event_start[id_in_read - self.window_size // 2]
      end = read.event_start[id_in_read + self.window_size // 2] + read.event_length[id_in_read + self.window_size // 2]
      signals = read.normalized_signal[start:end]
      self.update_base(reference, base, signals)