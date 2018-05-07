import h5py
from alphabet import *
import numpy as np
import numpy.ctypeslib
import copy
import statistics
import ctypes
import c_wrapper

def kmer_to_id(kmer):
  result = 0
  for c in kmer:
    result *= len(alphabet)
    result += inv_alphabet[c]
  return result


class tombo_kmer_model:
  def __init__(self, filename = None):
    self.k = None
    self.central_position = None
    self.mean = None
    self.sigma = None
    self.add = None
    if filename != None:
      self.load_from_hdf5(filename)

  def load_from_hdf5(self, filename):
    with h5py.File(filename, "r") as f:
      self.central_position = f.attrs["central_pos"]
      model = f["model"]
      kmer_model = {}
      for l in model:
        kmer_model[l[0].decode("ascii")] = [l[1], l[2]]
        self.k = len(l[0].decode("ascii"))
      self.mean = np.zeros(len(kmer_model))
      self.sigma = np.zeros(len(kmer_model))
      self.add = np.zeros(len(kmer_model))
      for kmer, value in kmer_model.items():
        id = kmer_to_id(kmer)
        self.mean[id] = value[0]
        self.sigma[id] = value[1]
        self.add[id] = np.log(1 / np.sqrt(2 * np.pi * self.sigma[id] ** 2))

  def median_filter(self, a, radius=7):
    res = copy.deepcopy(a)
    for i in range(radius // 2, len(a) - radius // 2):
      res[i] = statistics.median(a[i - radius // 2: i + radius // 2 + 1])
    return res

  def normalize_signal(self, raw_signal, file):
    normalization_meta = file['Analyses/RawGenomeCorrected_000/BaseCalled_template'].attrs
    scale = normalization_meta['scale']
    shift = normalization_meta['shift']
    return np.clip(self.median_filter((raw_signal - shift) / scale, 1), -4, 4)

  def log_chance_of_signal(self, kmer, value, stdv=None):
    if not isinstance(kmer, int):
      kmer = kmer_to_id(kmer)
    return self.add[kmer] - ((value - self.mean[kmer]) ** 2) / ((2 * self.sigma[kmer] ** 2))

  def get_c_object(self):
    return c_wrapper.new_KmerModel(self.k, self.central_position, self.mean, self.sigma)

  def free_c_object(self, obj):
    c_wrapper.delete_KmerModel(obj)


class picoamp_kmer_model:
  def __init__(self, filename=None):
    self.k = None
    self.central_position = None
    self.mean = None
    self.sigma = None
    self.sd_mean = None
    self.sd_lambda = None
    self.add = None
    if filename != None:
      self.load_from_file(filename)

  def rescale_value_to_tombo(self, x):
    return  x / 10.08832601 - 9.02667786411

  def rescale_difference_to_tombo(self, x):
    return x / 10.08832601

  def load_from_file(self, filename):
    with open(filename, "r") as f:
      kmer_map = {}
      for l in f:
        tokens = l.split()
        if tokens[0] == "kmer":
          continue
        kmer_map[tokens[0]] = list(map(float,tokens[1:]))
        self.k = len(tokens[0])
      self.central_position = self.k // 2
      self.mean = [None for i in kmer_map]
      self.sigma = [None for i in kmer_map]
      self.sd_mean = [None for i in kmer_map]
      self.sd_lambda = [None for i in kmer_map]
      self.add = [None for i in kmer_map]
      for kmer, value in kmer_map.items():
        id = kmer_to_id(kmer)
        self.mean[id] = self.rescale_value_to_tombo(value[0])
        self.sigma[id] = self.rescale_difference_to_tombo(value[1])
        self.sd_mean[id] = self.rescale_difference_to_tombo(value[2])
        self.sd_lambda[id] = (self.sd_mean[id] ** 3) / (self.rescale_difference_to_tombo(value[3]) ** 2)
        self.add[id] = np.log(1/np.sqrt(2 * np.pi * self.sigma[id]**2))

  def normalize_signal(self, raw_signal, file):
    normalization_meta = file['Analyses/RawGenomeCorrected_000/BaseCalled_template'].attrs
    scale = normalization_meta['scale']
    shift = normalization_meta['shift']
    return (raw_signal - shift) / scale

  #def normalize_signal(self, raw_signal, file):
  #  meta = file['UniqueGlobalKey/channel_id'].attrs
  #  return (raw_signal + meta["offset"]) * meta["range"] / meta["digitisation"]

  def log_chance_of_signal(self, kmer, value, stdv=None):
    if not isinstance(kmer, int):
      kmer = kmer_to_id(kmer)
    return self.add[kmer] - ((value - self.mean[kmer]) ** 2) / ((2 * self.sigma[kmer] ** 2))
