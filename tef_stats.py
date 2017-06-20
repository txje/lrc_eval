import sys
import argparse
import fasta

def main(unc, cor, fa, sorted=False, verbose=False):

  if verbose:
    print "Reading pacbio fasta"
  reads, names = fasta.read_fasta(fa)

  cor_aligned = 0
  unc_aligned = 0

  tp = 0
  fp = 0
  fn = 0
  #ne = 0
  tn = 0

  '''
  From ec_toolkit compute-stats.py:

  errorStats['TP'] += len(errPreCorrect.difference(errPostCorrect))
  errorStats['FP'] += len(errPostCorrect.difference(errPreCorrect))
  errorStats['FN'] += len(errPreCorrect.intersection(errPostCorrect))
  errorStats['NE'] += getNumWrongBase(errPreCorrect,errPostCorrect)

  # apparently, NE is the number of bases changed, but still incorrect
  '''

  '''
  From Error Correction Toolkit paper:

  We use the following measures for each program:
  number of erroneous bases identified and
  successfully corrected (true positives, TP), correct
  bases wrongly identified as errors and changed
  (false positives, FP), and erroneous bases that were
  either uncorrected or falsely corrected (false negatives,
  FN). We report sensitivity and specificity for
  each program. Then, we combine these into the gain
  metric [21], defined by gain = (TP - FP) /
  (TP + FN), which is the percentage of errors
  removed from the data set by the error-correction
  program. A negative gain value indicates that more
  errors have been introduced due to false corrections,
  which is not captured by measures such as sensitivity
  and specificity.
  '''

  if not sorted:
    uncor = {}
    corr = {}

    if verbose:
      print "Reading uncorrected TEF"
    for line in open(unc):
      data = line.strip().split(' ')
      fields = [int(a) for a in data[1:]]
      assert fields[0] == (len(fields)-1) / 4, "Number of errors does not match list"
      uncor[data[0]] = [fields[i:i+4] for i in xrange(1, len(fields), 4)]

    if verbose:
      print "Reading corrected TEF"
    for line in open(cor):
      data = line.strip().split(' ')
      fields = [int(a) for a in data[1:]]
      assert fields[0] == (len(fields)-1) / 4, "Number of errors does not match list"
      corr[data[0]] = [fields[i:i+4] for i in xrange(1, len(fields), 4)]

    if verbose:
      print "Some uncorrected reads:"
      print uncor.keys()[:10]
      print
      print "Some corrected reads:"
      print corr.keys()[:10]

    for n in names:
      if not uncor.has_key(n) and not corr.has_key(n):
        continue
      if not uncor.has_key(n) and corr.has_key(n):
        cor_aligned += 1
        continue
      if uncor.has_key(n) and not corr.has_key(n):
        unc_aligned += 1
        continue
      cor_aligned += 1
      unc_aligned += 1
      un = uncor[n]
      co = corr[n]

      if verbose:
        print
        print n
        print "%i errors in uncorrected read" % len(un)
        print "%i errors in corrected read" % len(co)

      cpos = set([c[0] for c in co])
      upos = set([u[0] for u in un])

      read_tp = len(upos - cpos)
      read_fp = len(cpos - upos)
      read_fn = len(cpos & upos)
      read_tn = len(reads[n]) - read_tp - read_fp - read_fn

      # not obviously trivial to compute this using set operations
      #read_ne = 0

      tp += read_tp
      fp += read_fp
      fn += read_fn
      tn += read_tn

  else: # sorted TEF
    cor_in = open(cor)
    uncor_in = open(unc)

    cor_aligned += 1
    cor_line = cor_in.readline()
    cor_data = cor_line.strip().split(' ')
    #cor_fields = [int(a) for a in cor_data[1:]]

    unc_aligned += 1
    uncor_line = uncor_in.readline()
    uncor_data = uncor_line.strip().split(' ')
    #uncor_fields = [int(a) for a in uncor_data[1:]]

    while len(cor_line) > 0 and len(uncor_line) > 0:

      if cor_data[0] == uncor_data[0]:
        n = cor_data[0]
        #co = [cor_fields[i:i+4] for i in xrange(1, len(cor_fields), 4)]
        #un = [uncor_fields[i:i+4] for i in xrange(1, len(uncor_fields), 4)]
        co = [cor_data[i] for i in xrange(1, len(cor_data), 4)]
        un = [uncor_data[i] for i in xrange(1, len(uncor_data), 4)]

        if verbose:
          print
          print n
          print "%i errors in uncorrected read" % len(un)
          print "%i errors in corrected read" % len(co)

        #cpos = set([c[0] for c in co])
        #upos = set([u[0] for u in un])
        cpos = set(co)
        upos = set(un)

        read_tp = len(upos - cpos)
        read_fp = len(cpos - upos)
        read_fn = len(cpos & upos)
        read_tn = len(reads[n]) - read_tp - read_fp - read_fn

        # not obviously trivial to compute this using set operations
        #read_ne = 0

        tp += read_tp
        fp += read_fp
        fn += read_fn
        tn += read_tn

      if uncor_data[0] <= cor_data[0]:
        unc_aligned += 1
        uncor_line = uncor_in.readline()
        if len(uncor_line) > 0:
          uncor_data = uncor_line.strip().split(' ')
        #uncor_fields = [int(a) for a in uncor_data[1:]]

      else: #if cor_data[0] < uncor_data[0]:
        cor_aligned += 1
        cor_line = cor_in.readline()
        if len(cor_line) > 0:
          cor_data = cor_line.strip().split(' ')
        #cor_fields = [int(a) for a in cor_data[1:]]


    if cor_line is None:
      while uncor_line is not None:
        unc_aligned += 1
        uncor_line = uncor_in.readline()
    if uncor_line is None:
      while cor_line is not None:
        cor_aligned += 1
        cor_line = cor_in.readline()
    cor_in.close()
    uncor_in.close()

  if tp + fn == 0:
    raise Exception("No read names matched (uncor: {}, cor: {})".format(uncor_data[0], cor_data[0]))

  gain = float(tp - fp) / (tp + fn)

  print "Sample\tMethod\tUncorrected reads\tCorrected reads\tRead gain/loss\tTP\tFP\tTN\tFN\tSensitivity\tSpecificity\tGain"
  print "%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%.4f\t%.4f\t%.4f" % (fa, cor, unc_aligned, cor_aligned, (cor_aligned-unc_aligned), tp, fp, tn, fn, (float(tp)/(tp + fn)), (float(tn)/(tn + fp)), gain)

  '''
  print "Uncorrected aligned reads: %i" % unc_aligned
  print "Corrected aligned reads: %i" % cor_aligned
  print "Change (positive is better): %i" % (cor_aligned - unc_aligned)
  print "TP: %i" % tp
  print "TN: %i" % tn
  print "FP: %i" % fp
  print "FN: %i" % fn
  #print "NE: %i" % ne
  print "Sensitivity (true positive rate): %.5f" % (float(tp)/(tp + fn))
  print "Specificity (true negative rate): %.5f" % (float(tn)/(tn + fp))
  print "Gain: %.2f" % gain
  '''


if __name__ == "__main__":
  parser = argparse.ArgumentParser("Compute true/false positive and gain stats for long-read error corrections")
  parser.add_argument("uncorrected", help="TEF file for original reads")
  parser.add_argument("corrected", help="TEF file for corrected reads")
  parser.add_argument("fasta", help="Uncorrected fasta")
  parser.add_argument("--verbose", help="Verbosity (every read)", default=False, action="store_true")
  parser.add_argument("--sorted", help="Allow sorted TEF to be read on the fly (not in memory)", default=False, action="store_true")
  args = parser.parse_args()
  main(args.uncorrected, args.corrected, args.fasta, args.sorted, args.verbose)
