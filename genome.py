from alphabet import complement

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
        current = genome(l.rstrip())
        result.append(current)
      else:
        current.append_line(l.rstrip())
  return result

def load_fastq(filename):
  with open(filename, 'r') as f:
    return parse_fastq(f)

def parse_fastq(lines):
  result = []
  state = 'balast'
  current = None
  for l in lines:
    if len(l) > 0 and l[0] == '@':
      current = genome(l)
      result.append(current)
      state = 'sequence'
    elif state == 'sequence':
      current.append_line(l.rstrip())
      state = 'balast'
  return result

def reverse_complement(r):
  res = [complement[x] for x in r]
  res.reverse()
  return res
