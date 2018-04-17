from read import resquiggled_read
from genome import parse_fastq

class basecalled_read(resquiggled_read):
  def __init__(self, filename = None, kmer_model = None, basecall_group = None):
    self.basecall_group = basecall_group
    self.fastq = None
    resquiggled_read.__init__(self, filename, kmer_model)

  def load_from_open_fast5(self, file, kmer_model):
    resquiggled_read.load_from_open_fast5(self, file, kmer_model)
    if self.basecall_group != None:
      fastq = file[self.basecall_group + "/BaseCalled_template/Fastq"].value.decode('ascii')
      self.basecall = parse_fastq(fastq.split('\n'))[0].bases

