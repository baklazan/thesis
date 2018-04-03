import h5py
import os
import shutil
import argparse
from genome import *

parser = argparse.ArgumentParser()
parser.add_argument("read_basedir", help="base directory of resquiggled fast5 files")
parser.add_argument("output_dir", help="directory to copy files to (must exist)")
parser.add_argument("reference", help="fasta file with reference genome")
parser.add_argument("begin", help="begining of interval (inclusive)")
parser.add_argument("end", help="end of interval (exclusive)")


args = parser.parse_args()

if args.read_basedir == args.output_dir:
  print("read_basedir and output_dir must be different")
  exit(1)

begin = int(args.begin)
end = int(args.end)
if end - begin <= 0:
  print("interval must have positive length")
  exit(1)

read_files = [file for file in os.listdir(args.read_basedir) if not os.path.isdir(os.path.join(args.read_basedir, file))]
read_files = filter(lambda x : x[-6:] == ".fast5", read_files)

reference = load_fasta(args.reference)[0].bases

for read_file in read_files:
  filename = os.path.join(args.read_basedir, read_file)
  good = False
  print(filename)
  try:
    with h5py.File(filename, "r") as f:
      read_id = list(f['Raw/Reads'].keys())[0]

      alignment_meta = f['Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment'].attrs
      start_in_reference = alignment_meta['mapped_start']
      end_in_reference = alignment_meta['mapped_end']
      strand = alignment_meta['mapped_strand'].decode('ascii')
      if strand == '-':
        start_in_reference, end_in_reference = len(reference) - end_in_reference, len(reference) - start_in_reference
      if min(end_in_reference, end) - max(start_in_reference, begin) > 0:
        good = True

  except KeyError:
    continue
  if good:
    copy_name = os.path.join(args.output_dir, read_file)
    print("copying {} to {}".format(filename, copy_name))
    shutil.copy2(filename, copy_name)
