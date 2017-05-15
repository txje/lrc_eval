import argparse
import aln_formats

def main(sam, tef, rename, fa):

  read_names = None
  if rename in ["sprai"]:
    if fa is None:
      raise Exception("A matching uncorrected fasta file is required to convert sprai indexed names back to their original")
    import fasta
    _, read_names = fasta.read_fasta(fa)

  fout = open(tef, 'w')
  l = 0

  for query_name, alignments in aln_formats.iter_maf(maf, 0, 0, 1000000000, by_query=True):
    al = sorted(alignments, key = lambda al: (al.accuracy() * abs(al.query.end - al.query.start)))[-1] # keep only "best" alignment, by total correct bp
    #fout.write('\n' + str(al))

    '''
TEF (format)

readid n-errors [pos tb wb ind]+

In the above format, the fields are described as below :

Fields    Description
readid    ID of the read corrected
n-errors  Integer. Number of errors corrected in the read.
pos       Position for fix (0 < = pos < length of the read)
tb        true value of the base at pos.
wb        wrong value of the base at pos.
          wb should be current base at read
          tb,wb is one of {0,1,2,3,4,5}
          0 = 'A', 1 = 'C', 2 = 'G', 3 = 'T', 5 = '-'
ind       indicates the type of error. one of {0,1,2}
          0 substitution (bad char in the read at pos) or
          1 deletion (missing char in the read after pos) or
          2 insertion (extra char in the read at pos)
    '''

    n_errs = 0
    err_strings = []
    q = al.query.start - 1
    t = al.target.start - 1
    # whole bunch of ambiguity codes will all map to N (they are present in the a_thaliana reference...)
    almap = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4, '-':5, 'R':4, 'Y':4, 'S':4, 'W':4, 'K':4, 'M':4, 'B':4, 'V':4, 'D':4, 'H':4}
    for i in xrange(len(al.query.alignment)):
      qlocus = al.query.alignment[i]
      tlocus = al.target.alignment[i]
      if qlocus == tlocus:
        continue

      n_errs += 1
      if tlocus == '-':
        ind = 2
        q += 1
      elif qlocus == '-':
        ind = 1
        t += 1
      else:
        ind = 0
        q += 1
        t += 1
      err_strings.append("%i %i %i %i" % (q, almap[tlocus.upper()], almap[qlocus.upper()], ind))

    if rename is not None:
      if rename in ["ectools"]:
        query_name = query_name[:query_name.rindex("_corrected")]
      elif rename == "sprai":
        query_name = read_names[int(query_name[:query_name.index('/')])]
      elif rename in ["nanocorr"]:
        query_name = query_name[:query_name.rindex("_consensus")]

    tef_line = "%s %i %s" % (query_name, n_errs, ' '.join(err_strings))
    fout.write(('\n' if l > 0 else '') + tef_line)
    l += 1

  fout.close()


if __name__ == "__main__":
  parser = argparse.ArgumentParser("Convert SAM to TEF format (I hate bioinformatics)")
  parser.add_argument("sam", help="SAM alignment file")
  parser.add_argument("tef", help="Output TEF file")
  parser.add_argument("--rename", help="Method used to rename reads, if any")
  parser.add_argument("--fasta", help="Original uncorrected FASTA file, for names mostly")
  args = parser.parse_args()
  main(args.sam, args.tef, args.rename, args.fasta)
