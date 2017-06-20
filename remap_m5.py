import fasta
import argparse

def main(m5, fa, method, trim_suffix=False):
  read_names = None
  if method in ["sprai"]:
    if fa is None:
      raise Exception("A matching uncorrected fasta file is required to convert sprai indexed names back to their original")
    import fasta
    _, read_names = fasta.read_fasta(fa)

  for line in open(m5):
    fields = line.strip().split()
    n = fields[0]

    if '/' in n and trim_suffix:
      n = n[:n.rindex('/')]

    if method is not None:
      if method in ["ectools"]:
        n = n[:n.rindex("_corrected")]
      elif method == "sprai":
        n = read_names[int(n[:n.index('/')]) - 1] # !! reads are 1-indexed as of sprai v0.9.9.23
      elif method in ["nanocorr"]:
        n = n[:n.rindex("_consensus")]
      elif method in ["lorma"]:
        if '_' in n:
          n = n[:n.rindex("_")] # lorma adds a "_<index>" to each read, which increments from 1 for each subread it was split into
    fields[0] = n
    print ' '.join(fields)

if __name__ == "__main__":
  parser = argparse.ArgumentParser("Map names in m5 file back to original (uncorrected) names")
  parser.add_argument("m5", help="m5 file for corrected reads")
  parser.add_argument("fasta", help="Uncorrected fasta")
  parser.add_argument("method", help="Method used to correct reads - each changes the names in its own way")
  parser.add_argument("--trim_suffix", help="The input m5 file has the blasr '/start_end' suffix", default=False, action="store_true")
  args = parser.parse_args()
  main(args.m5, args.fasta, args.method, args.trim_suffix)
