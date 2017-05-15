
from string import maketrans
COMPL = maketrans("ATGC","TACG")

def rc(seq):
  return seq.translate(COMPL)[::-1]

def read_fasta(fasta):
  data = open(fasta).read().strip().split('\n')
  reads = {}
  names = []
  name = None
  seq = ""
  for d in data:
    if len(d) == 0:
      continue
    if d[0] == '>':
      if name is not None:
        reads[name] = seq
      name = d[1:]
      if ' ' in name:
        name = name[:name.index(' ')]
      names.append(name)
      seq = ""
    else:
      seq += d
  reads[name] = seq
  return reads, names

def write_fasta(fname, reads, names=None):
  fout = open(fname, 'w')
  if names is None:
    names = reads.keys()
  for r in xrange(len(names)):
    if r > 0:
      fout.write('\n')
    fout.write(">%s" % names[r])
    seq = reads[names[r]]
    fout.write("\n%s" % '\n'.join([seq[i:min(len(seq),i+100)] for i in xrange(0, len(seq), 100)]))
  fout.close()

