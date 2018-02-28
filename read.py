import h5py
import numpy as np
from kmer_model import *
import matplotlib.pyplot as plt
from scipy import interpolate


class resquiggled_read:
  def __init__(self, filename = None, kmer_model = None):
    self.raw_signal = None
    self.normalized_signal = None
    self.start_in_reference = None
    self.end_in_reference = None
    self.event_start = None
    self.event_length = None
    self.event_mean = None
    self.event_base = None
    self.number_of_events = None
    self.read_id = None
    if filename != None:
      self.load_from_fast5(filename, kmer_model)
  
  def load_from_fast5(self, filename, kmer_model):
    with h5py.File(filename, 'r') as f:
      self.read_id = list(f['Raw/Reads'].keys())[0]
      self.raw_signal = f['Raw/Reads/{}/Signal'.format(self.read_id)].value
      
      alignment_meta = f['Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment'].attrs
      self.start_in_reference = alignment_meta['mapped_start']
      self.end_in_reference = alignment_meta['mapped_end']
      
      
      events = f['Analyses/RawGenomeCorrected_000/BaseCalled_template/Events']
      relative_start = events.attrs['read_start_rel_to_raw']
      events = events.value
      self.event_mean = np.array([e[0] for e in events])
      self.event_start = np.array([e[2] + relative_start for e in events])
      self.event_length = np.array([e[3] for e in events])
      self.event_base = np.array([e[4] for e in events])
      
      self.normalized_signal = kmer_model.normalize_signal(self.raw_signal, f)
      self.number_of_events = len(events)

  def tweak_normalization(self, reference, kmer_model):
    data = []
    for i, mean in enumerate(self.event_mean):
      id_in_reference = i + self.start_in_reference
      kmer_start = id_in_reference - kmer_model.central_position
      kmer_end = kmer_start + kmer_model.k
      kmer = reference[kmer_start:kmer_end]
      expected_mean = kmer_model.mean[kmer_to_id(kmer)]
      if abs(expected_mean - mean) <= 1:
        data.append((mean, expected_mean))

    data.sort()
    x = [d[0] for d in data]
    y = [d[1] for d in data]
    fnc = interpolate.splrep(x, y, s=len(x))
    self.normalized_signal = interpolate.splev(self.normalized_signal, fnc)
