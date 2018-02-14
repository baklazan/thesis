import h5py
from alphabet import *
import numpy as np

class kmer_model:
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
      self.mean = [None for i in kmer_model]
      self.sigma = [None for i in kmer_model]
      self.add = [None for i in kmer_model]
      for kmer, value in kmer_model.items():
        id = self.kmer_to_id(kmer)
        self.mean[id] = value[0]
        self.sigma[id] = value[1]
        self.add[id] = np.log(1/np.sqrt(2 * np.pi * value[1]**2))
    
  def kmer_to_id(self, kmer):
    result = 0
    for c in kmer:
      result *= len(alphabet)
      result += inv_alphabet[c]
    return result
    
  def log_chance_of_signal(self, kmer, value):
    if not isinstance(kmer, int):
      kmer = self.kmer_to_id(kmer)
    return self.add[kmer] - ((value-self.mean[kmer])**2)/((2 * self.sigma[kmer]**2))