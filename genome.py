class genome:
  def __init__(self, desc_line = None):
    self.bases = []
    self.description = desc_line
  def append_line(self, l):
    for c in l:
      self.bases.append(c)

def load_fasta(filename):
  result = []
  with open(filename, 'r') as f:
    current = None
    for l in f:
      if l[0] == '>':
        current = genome(l)
        result.append(current)
      else:
        current.append_line(l.rstrip())
  return result
  