from kmer_model import *
from base import *
import random
import matplotlib.pyplot as plt

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

def log_average(a, b):
  if np.isinf(a) and np.isinf(b):
    return -np.isinf
  if a < b:
    return b + np.log(1 + np.exp(a - b)) - np.log(2)
  else:
    return a + np.log(1 + np.exp(b - a)) - np.log(2)

class window_model:
  def __init__(self, kmer_model = None, window_size = 13, min_event_length = 3, op = None, buffer_size = 1, penalty = 1.0):
    self.window_size = window_size
    self.min_event_length = min_event_length
    self.kmer_model = kmer_model
    self.op = op
    if op == None:
      self.op = max_operation()
    self.buffer_size = buffer_size
    self.penalty = penalty

  def reresquiggle(self, signal, reference_kmer_ids, start_buffer=0, end_buffer=0, plot = None):
    dp = np.full((len(signal) + 1, len(reference_kmer_ids)*2), self.op.neutral_element(), dtype=float)
    reconstruct_length = np.full((len(signal) + 1, len(reference_kmer_ids)*2), -1, dtype=int)
    for j in range(start_buffer+1):
      dp[j][0] = self.op.start_element()-j*self.penalty
    for i, kmer_id in enumerate(reference_kmer_ids):
      for j in range(self.min_event_length, len(signal) + 1):
        accumulator = 0
        for l in range(1, self.min_event_length + 1):
          p = self.kmer_model.log_chance_of_signal(kmer_id, signal[j-l])
          accumulator = self.op.append_combine(p, accumulator)
        val_if_shortest = self.op.append_combine(dp[j-self.min_event_length][i*2], accumulator)
        if val_if_shortest > dp[j][i*2 + 1]:
          reconstruct_length[j][i*2+1] = self.min_event_length
        dp[j][i*2 + 1] = self.op.combine(dp[j][i*2+1], val_if_shortest)
        p = self.kmer_model.log_chance_of_signal(kmer_id, signal[j-1])
        val_if_long = self.op.append_combine(dp[j-1][i*2+1], p)
        if val_if_long > dp[j][i*2+1]:
          reconstruct_length[j][i*2+1] = reconstruct_length[j-1][i*2+1] + 1
        dp[j][i*2 + 1] = self.op.combine(dp[j][i*2+1], val_if_long)
      if i+1 < len(reference_kmer_ids):
        for j in range(len(signal) + 1):
          val_if_no_flashback = dp[j][i*2+1]
          if val_if_no_flashback > dp[j][i*2+2]:
            reconstruct_length[j][i*2+2] = 0
          dp[j][i*2+2] = self.op.combine(dp[j][i*2+2], val_if_no_flashback)
          if j >= 1:
            p1 = self.kmer_model.log_chance_of_signal(kmer_id, signal[j-1])
            p2 = self.kmer_model.log_chance_of_signal(reference_kmer_ids[i+1], signal[j-1])
            p = log_average(p1, p2)
            val_if_flashback = self.op.append_combine(dp[j-1][i*2+2], p)
            if val_if_flashback > dp[j][i*2+2]:
              reconstruct_length[j][i*2+2] = reconstruct_length[j-1][i*2+2] + 1
            dp[j][i*2+2] = self.op.combine(dp[j][i*2+2], val_if_flashback)

    for i in range(end_buffer+1):
      dp[-i-1][-1] -= i*self.penalty
    res = np.max(dp[-end_buffer-1:, -1])
    if plot != None:
      events = []
      cur = np.argmax(dp[-end_buffer-1:, -1]) + len(signal) -end_buffer
      for i in range(len(reference_kmer_ids)*2-1, 0, -1):
        l = reconstruct_length[cur][i]
        events.append((cur - l, l))
        cur -= l
      x = []
      y = []
      events.reverse()
      for i, e in enumerate(events[0::2]):
        for j in range(e[0], e[0]+e[1]):
          x.append(j)
          y.append(self.kmer_model.mean[reference_kmer_ids[i]])

      plot.plot(signal)
      plot.plot(x, y)
    return res

  def update_base(self, reference, base, signal_values, start_buffer, end_buffer):
    ref_st = base.id - self.window_size // 2 - self.kmer_model.central_position
    ref_en = ref_st + self.window_size + self.kmer_model.k - 1
    if ref_st < 0 or ref_en > len(reference):
      return
    context = list(reference[ref_st:ref_en])
    for i, c in enumerate(alphabet):
      context[self.kmer_model.central_position + self.window_size // 2] = c
      kmer_ids = [kmer_to_id(context[j:j + self.kmer_model.k]) for j in range(len(context) - self.kmer_model.k + 1)]
      base.log_probability[i] += self.reresquiggle(signal_values, kmer_ids, start_buffer, end_buffer)

  def show_base(self, reference, read, base):
    id_in_read = base.id - read.start_in_reference
    kmer_start_id = id_in_read - self.window_size // 2
    kmer_end_id = id_in_read + self.window_size // 2
    signal_start_id = kmer_start_id - self.buffer_size
    signal_end_id = kmer_end_id + self.buffer_size
    if signal_start_id < 0 or signal_end_id >= len(read.event_start):
      return
    signal_start = read.event_start[signal_start_id]
    signal_end = read.event_start[signal_end_id] + read.event_length[signal_end_id]
    signals = read.normalized_signal[signal_start:signal_end]
    start_buffer = read.event_start[signal_start_id + 2 * self.buffer_size] - signal_start
    end_buffer = signal_end - (read.event_start[signal_end_id - 2 * self.buffer_size] + read.event_length[
      signal_end_id - 2 * self.buffer_size])
    context_st = base.id - self.window_size // 2 - self.kmer_model.central_position
    context_en = context_st + self.window_size + self.kmer_model.k - 1
    if context_st < 0 or context_en > len(reference):
      return
    context = list(reference[context_st:context_en])
    rows = np.floor(np.sqrt(len(alphabet)))
    columns = np.ceil(len(alphabet) / rows)
    plt.figure()
    plt.suptitle("{}: {}<-{}".format(base.id, base.reference_value, base.real_value))
    plots = [plt.subplot(rows, columns, i + 1) for i in range(len(alphabet))]
    for i, c in enumerate(alphabet):
      context[self.kmer_model.central_position + self.window_size // 2] = c
      kmer_ids = [kmer_to_id(context[j:j + self.kmer_model.k]) for j in range(len(context) - self.kmer_model.k + 1)]
      plots[i].set_title("{} {}".format(c, base.log_probability[i]))
      self.reresquiggle(signals, kmer_ids, start_buffer, end_buffer, plots[i])
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()


  def update_probabilities(self, reference, read, interesting_bases):
    for base in interesting_bases:
      id_in_read = base.id - read.start_in_reference
      kmer_start_id = id_in_read - self.window_size // 2
      kmer_end_id = id_in_read + self.window_size // 2
      signal_start_id = kmer_start_id - self.buffer_size
      signal_end_id = kmer_end_id + self.buffer_size
      if signal_start_id < 0 or signal_end_id >= len(read.event_start):
        continue
      signal_start = read.event_start[signal_start_id]
      signal_end = read.event_start[signal_end_id] + read.event_length[signal_end_id]
      signals = read.normalized_signal[signal_start:signal_end]
      start_buffer = read.event_start[signal_start_id + 2*self.buffer_size] - signal_start
      end_buffer = signal_end - (read.event_start[signal_end_id - 2*self.buffer_size] + read.event_length[signal_end_id - 2*self.buffer_size])
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
