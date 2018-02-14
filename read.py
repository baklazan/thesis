import h5py
import numpy as np

class resquiggled_read:
  def __init__(self, filename = None):
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
      self.load_from_fast5(filename)
  
  def load_from_fast5(self, filename):
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
      
      normalization_meta = f['Analyses/RawGenomeCorrected_000/BaseCalled_template'].attrs
      scale = normalization_meta['scale']
      shift = normalization_meta['shift']
      self.normalized_signal = (self.raw_signal - shift) / scale
      
      self.number_of_events = len(events)
