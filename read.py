import h5py
import numpy as np
from kmer_model import *
import matplotlib.pyplot as plt
from scipy import interpolate
import c_wrapper


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
    self.strand = None
    if filename != None:
      self.load_from_fast5(filename, kmer_model)
  
  def load_from_fast5(self, filename, kmer_model):
    with h5py.File(filename, 'r') as f:
      self.load_from_open_fast5(f, kmer_model)

  def load_from_open_fast5(self, file, kmer_model):
    self.read_id = list(file['Raw/Reads'].keys())[0]
    self.raw_signal = file['Raw/Reads/{}/Signal'.format(self.read_id)].value

    alignment_meta = file['Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment'].attrs
    self.start_in_reference = alignment_meta['mapped_start']
    self.end_in_reference = alignment_meta['mapped_end']
    self.strand = alignment_meta['mapped_strand'].decode('ascii')

    events = file['Analyses/RawGenomeCorrected_000/BaseCalled_template/Events']
    relative_start = events.attrs['read_start_rel_to_raw']
    events = events.value
    self.event_mean = np.array([e[0] for e in events])
    self.event_start = np.array([e[2] + relative_start for e in events], dtype=np.int32)
    self.event_length = np.array([e[3] for e in events], dtype=np.int32)
    self.event_base = np.array([e[4] for e in events])

    self.normalized_signal = kmer_model.normalize_signal(self.raw_signal, file)
    self.number_of_events = len(events)
    self.already_fixed = False

  def fix_start_in_reference(self, reference):
    if self.strand == '-' and not self.already_fixed:
      self.start_in_reference, self.end_in_reference = len(reference) - self.end_in_reference, len(reference) - self.start_in_reference
    self.already_fixed = True

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

  def get_c_object(self):
    return c_wrapper.new_ResquiggledRead(self.normalized_signal,
                                        len(self.normalized_signal),
                                        self.event_start,
                                        self.event_length,
                                        self.number_of_events,
                                        self.start_in_reference)

  def free_c_object(self, obj):
    c_wrapper.delete_ResquiggledRead(obj)
