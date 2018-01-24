import h5py
import numpy as np
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", help="output file for a posteriori probabilities")
parser.add_argument("-i", "--interesting", help="list of interesting positions in reference")
parser.add_argument("-p", "--print_best", help="print top n candidates")
parser.add_argument("-k", "--kmer_model", help="file with kmer model to use")
parser.add_argument("reference", help="reference fasta file")
parser.add_argument("read", help="resquiggled fast5 file(s)", nargs="+")
args = parser.parse_args()

alphabet = ['A', 'C', 'G', 'T']
inv_alphabet = {c : i for i, c in enumerate(alphabet)}
window_size = 11
min_event_length = 2
max_event_length = 30

def kmer_to_int(kmer):
  result = 0
  for c in kmer:
    result *= len(alphabet)
    result += inv_alphabet[c]
  return result

def convert_to_list(model):
  res = [None for i in range(len(model))]
  for key, val in model.items():
    res[kmer_to_int(key)] = val + [np.log(1/np.sqrt(2 * np.pi * val[1]**2))]
  return res

def normal_distribution(x, mean, sigma):
  return np.exp(-((x - mean)**2)/(2*sigma**2)) / np.sqrt(2*np.pi*sigma**2)

def log_chance_of_signal(kmer_id, val):
  mean, sigma, add = kmer_model[kmer_id][:3]
  return add - ((val-mean)**2)/((2 * sigma**2))
  #return np.log(normal_distribution(val, mean, sigma))

def reresquiggle(signal, reference_kmer_ids):
  prefix_sums = np.zeros(len(signal)+1)
  for i, s in enumerate(signal):
    prefix_sums[i+1] = prefix_sums[i] + signal[i]
  dp = np.full((len(signal)+1, len(reference_kmer_ids)+1), -np.inf, dtype = float)
  dp[0][0] = 0
  for i, kmer_id in enumerate(reference_kmer_ids):
    for j in range(len(signal)+1):
      for l in range(min_event_length, min(max_event_length, j)+1):
        average = (prefix_sums[j] - prefix_sums[j-l])/l
        dp[j][i+1] = max(dp[j][i+1], log_chance_of_signal(kmer_id, average) + dp[j-l][i])
  return dp[len(signal)][len(reference_kmer_ids)]

def update_probabilities(reference, index, signal_values, probabilities, dirty):
  ref_st = index-window_size//2-central_position
  ref_en = ref_st + window_size + k - 1
  context = list(reference[ref_st:ref_en])
  for i, c in enumerate(alphabet):
    context[central_position + window_size//2] = c
    kmer_ids = [kmer_to_int(context[j:j+k]) for j in range(len(context) - k - 1)]
    probabilities[index, i] = reresquiggle(signal_values, kmer_ids)
  dirty[index] = True



def process_fast5(filename, reference, probabilities, dirty):
  with h5py.File(filename, "r") as f:
    
    events = f['Analyses/RawGenomeCorrected_000/BaseCalled_template/Events']
    alignment_meta = f['Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment'].attrs
    mapped_start = alignment_meta['mapped_start']
    mapped_end = alignment_meta['mapped_end']
    print("[{}, {})".format(mapped_start, mapped_end))
    
    normalization_meta = f['Analyses/RawGenomeCorrected_000/BaseCalled_template'].attrs
    scale = normalization_meta['scale']
    shift = normalization_meta['shift']
    
    relative_start = f['Analyses/RawGenomeCorrected_000/BaseCalled_template/Events'].attrs['read_start_rel_to_raw']
    meta = f["UniqueGlobalKey/channel_id"].attrs
    read_id = list(f['Raw/Reads'].keys())[0]
    raw = f['Raw/Reads/{}/Signal'.format(read_id)]
    normalized_raw = (raw - shift) / scale
    
    for i in range(window_size//2, len(events)-window_size//2):
      start = events[i-window_size//2][2] + relative_start
      end = events[i+window_size//2][2] + events[i+window_size//2][3] + relative_start
      signals = normalized_raw[start:end]
      update_probabilities(reference, mapped_start + i, signals, probabilities, dirty)
      if(i % 100 == 0):
        print("{} / {}".format(i, len(events)-k+1))

def normalize(probs, dirty = None):
  indices = range(len(probs))
  if not dirty is None:
    indices = [i for i, v in enumerate(dirty) if v]
  for i in indices:
    sum = 0
    for j, v in enumerate(probs[i]):
      probs[i][j] = np.exp(v)
      sum += probs[i][j]
    for j, v in enumerate(probs[i]):
      probs[i][j] /= sum

def print_entry(id, bullet, reference, probabilities, interesting):
  if id not in interesting:
    bullet = " "
  sys.stdout.write("{} {}:\t{}".format(bullet, id, reference[id]))
  if id in interesting:
    sys.stdout.write("<-{}".format(interesting[id]))
  for i, p in enumerate(probabilities[id]):
    decoration = " "
    if alphabet[i] == reference[id]:
      decoration = "*"
    elif id in interesting and interesting[id] == alphabet[i]:
      decoration = "-"
    sys.stdout.write("\t{:10.8f}{}".format(p, decoration))
  print()

kmer_model_file = "data/tombo.DNA.model"
if args.kmer_model:
  kmer_model_file = args.kmer_model

k = None
central_position = None
kmer_model = {}
with h5py.File(kmer_model_file, "r") as f:
  model = f["model"]
  for l in model:
    kmer_model[l[0].decode("ascii")] = [l[1], l[2]]
    k = len(l[0].decode("ascii"))
  central_position = f.attrs["central_pos"]
kmer_model = convert_to_list(kmer_model)

reference = ""
reference_file = args.reference
with open(reference_file, "r") as f:
  lines = []
  first_started = False
  for l in f:
    if l[0] == '>':
      if first_started:
        break
      else:
        first_started = True
    else:
      lines.append(l.rstrip())
  reference = "".join(lines)

reference_kmer_id = np.zeros(len(reference))
cur = 0
number_of_kmers = len(alphabet)**k
for i, c in enumerate(reference):
  cur *= len(alphabet)
  cur += inv_alphabet[c]
  cur %= number_of_kmers
  if i >= k-1:
    reference_kmer_id[i -k - 1 + central_position] = cur

probabilities = np.zeros((len(reference), len(alphabet)))
dirty = np.full(len(reference), False, dtype=bool)
for read_file in args.read:
  process_fast5(read_file, reference, probabilities, dirty)

normalize(probabilities, dirty)

if args.output:
  with open(args.output, "w") as f:
    for i, p in enumerate(probabilities):
      f.write("{} {}\n".format(reference[i], "\t".join(map(str, p))))

interesting = {}
if args.interesting:
  with open(args.interesting, "r") as f:
    for l in f:
      tokens = l.split()
      interesting[int(tokens[0])] = tokens[1]
  for id, c in interesting.items():
    print_entry(id, " ", reference, probabilities, interesting)
  print()

if args.print_best:
  by_SNP_prob = [(1 - p[inv_alphabet[reference[i]]], i) for i, p in enumerate(probabilities) if dirty[i]]
  by_SNP_prob.sort()
  by_SNP_prob.reverse()
  for p, id in by_SNP_prob[:int(args.print_best)]:
    print_entry(id, "*", reference, probabilities, interesting)
    

count_good = 0
count_bad = 0
for i, p in enumerate(probabilities):
  if dirty[i]:
    if p[inv_alphabet[reference[i]]] == np.max(p):
      count_good += 1
    else:
      count_bad += 1

print("good: {} / {} [{} %]".format(count_good, count_good+count_bad, count_good/(count_good + count_bad)))