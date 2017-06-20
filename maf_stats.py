import argparse
import aln_formats
import fasta

def main(unc, cor, fa):

  print "Reading pacbio fasta"
  reads, names = fasta.read_fasta(fa)

  uncor = {}
  corr = {}

  print "Reading uncorrected alignments"
  for query_name, alignments in aln_formats.iter_maf(unc, 0, 0, 1000000000, by_query=True):
    al = sorted(alignments, key = lambda al: (al.accuracy() * abs(al.query.end - al.query.start)))[-1] # keep only "best" alignment, by total correct bp
    uncor[al.query.name] = al

  print "Reading corrected alignments"
  for query_name, alignments in aln_formats.iter_maf(cor, 0, 0, 1000000000, by_query=True):
    al = sorted(alignments, key = lambda al: (al.accuracy() * abs(al.query.end - al.query.start)))[-1]
    corr[al.query.name] = al

  new_aligned = 0
  new_unaligned = 0

  for n in names:
    if not uncor.has_key(n) and not corr.has_key(n):
      continue
    if not uncor.has_key(n) and corr.has_key(n):
      new_aligned += 1
      continue
    if uncor.has_key(n) and not corr.has_key(n):
      new_unaligned += 1
      continue
    un = uncor[n]
    co = corr[n]
    print
    print n
    print un
    print co
    for i in xrange(len(un.query.alignment)):
      qlocus = un.query.alignment[i]
      tlocus = un.target.alignment[i]
      if qlocus == tlocus:
        continue

      n_errs += 1
      if tlocus == '-':
        q += 1
      elif qlocus == '-':
        t += 1
      else:
        q += 1
        t += 1


if __name__ == "__main__":
  parser = argparse.ArgumentParser("Compute true/false positive and gain stats for long-read error corrections")
  parser.add_argument("uncorrected", help="MAF alignment file (a la LAST) for original reads")
  parser.add_argument("corrected", help="MAF alignment file for corrected reads")
  parser.add_argument("fasta", help="Uncorrected fasta")
  args = parser.parse_args()
  main(args.uncorrected, args.corrected, args.fasta)
