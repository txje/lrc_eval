import sys
import argparse
import fasta
import aln_formats


def main(unc, cor, fa, method, verbose=False):

  if verbose:
    print "Reading pacbio fasta"
  pacbio_reads, names = fasta.read_fasta(fa)

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

  cor_iter = aln_formats.iter_m5(cor, 0, 0, 1000000000, by_query=True)
  uncor_iter = aln_formats.iter_m5(unc, 0, 0, 1000000000, by_query=True)

  cor_query = None
  uncor_query = None

  correct_uncorrected = 0
  incorrect_uncorrected = 0
  correct_corrected = 0
  incorrect_corrected = 0

  # ------ !! keep track of loci on the target sequence since coordinates on the query sequence will change dramatically ------
  # ------ !! although if sequences align to different places, everything will go to crap ------

  while True:
    if cor_query is not None and cor_query == uncor_query and cor_best_aln.target.name == uncor_best_aln.target.name:

      if verbose:
        print
        print "{} aligned for both uncorrected and corrected".format(cor_query)
        print "{} errors in uncorrected read".format(len(incorrect_loci_in_uncorrected))
        print "{} errors in corrected read".format(len(incorrect_loci_in_corrected))

      read_tp = len(prev_incorrect & now_correct)
      read_fp = len(prev_correct & now_incorrect)
      read_fn = len(prev_incorrect & now_incorrect)
      read_tn = len(prev_correct & now_correct)

      # not obviously trivial to compute this using set operations
      #read_ne = 0

      tp += read_tp
      fp += read_fp
      fn += read_fn
      tn += read_tn

    if (uncor_query is None or cor_iter is None or uncor_query <= cor_query) and uncor_iter is not None:
      try:
        uncor_query, uncor_aln = uncor_iter.next()
      except:
        uncor_iter = None
        if cor_iter is None:
          break
      unc_aligned += 1
      uncor_best_aln = sorted(uncor_aln, key = lambda al: (al.accuracy() * abs(al.query.end - al.query.start)))[-1] # keep only "best" alignment, by total correct bp
      incorrect_loci_in_uncorrected = []
      correct_loci_in_uncorrected = []

      # get all incorrect loci in uncorrected alignment
      tpos = uncor_best_aln.target.start
      for i in xrange(len(uncor_best_aln.alignment)):
        if uncor_best_aln.alignment[i] == '|':
          correct_loci_in_uncorrected.append(tpos)
        else:
          if uncor_best_aln.target.alignment[i] == '-':
            incorrect_loci_in_uncorrected.append("{}i{}".format(uncor_best_aln.query.alignment[i], tpos)) # <nucleotide> inserted before tpos
          else:
            incorrect_loci_in_uncorrected.append("x{}".format(tpos)) # <nucleotide> mismatch or deleted at tpos
        if uncor_best_aln.target.alignment[i] != '-':
          tpos += 1
      incorrect_loci_in_uncorrected.extend(range(uncor_best_aln.target.start-uncor_best_aln.query.start, uncor_best_aln.target.start) + range(uncor_best_aln.target.end, uncor_best_aln.target.end+uncor_best_aln.query.length-uncor_best_aln.query.end)) # finagle an estimate of the target regions that are supposed to be covered by the read

      prev_incorrect = set(incorrect_loci_in_uncorrected)
      prev_correct = set(correct_loci_in_uncorrected)
      incorrect_uncorrected += len(prev_incorrect)
      correct_uncorrected += len(prev_correct)

    else:
      try:
        cor_query, cor_aln = cor_iter.next()
      except:
        cor_iter = None
        if uncor_iter is None:
          break
      cor_aligned += 1
      cor_best_aln = sorted(cor_aln, key = lambda al: (al.accuracy() * abs(al.query.end - al.query.start)))[-1] # keep only "best" alignment, by total correct bp
      incorrect_loci_in_corrected = []
      correct_loci_in_corrected = []

      # get all incorrect loci in corrected alignment
      tpos = cor_best_aln.target.start
      for i in xrange(len(cor_best_aln.alignment)):
        if cor_best_aln.alignment[i] == '|':
          correct_loci_in_corrected.append(tpos)
        else:
          if cor_best_aln.target.alignment[i] == '-':
            incorrect_loci_in_corrected.append("{}i{}".format(cor_best_aln.query.alignment[i], tpos)) # <nucleotide> inserted before tpos
          else:
            incorrect_loci_in_corrected.append("{}".format(tpos)) # <nucleotide> mismatch or deleted at tpos
        if cor_best_aln.target.alignment[i] != '-':
          tpos += 1
      incorrect_loci_in_corrected.extend(range(cor_best_aln.target.start-cor_best_aln.query.start, cor_best_aln.target.start) + range(cor_best_aln.target.end, cor_best_aln.target.end+cor_best_aln.query.length-cor_best_aln.query.end)) # finagle an estimate of the target regions that are supposed to be covered by the read

      now_incorrect = set(incorrect_loci_in_corrected)
      now_correct = set(correct_loci_in_corrected) # these should already be unique, but we need them to be sets to do set operations
      incorrect_corrected += len(now_incorrect)
      correct_corrected += len(now_correct)


  if tp + fn == 0:
    raise Exception("No read names matched (uncor: {}, cor: {})".format(uncor_query, cor_query))

  gain = float(tp - fp) / (tp + fn)

  print "Sample\tMethod\tUncorrected reads\tCorrected reads\tRead gain/loss\tUncorrected wrong bp\tUncorrected right bp\tCorrected wrong bp\tCorrected right bp\tTP\tFP\tTN\tFN\tSensitivity\tSpecificity\tGain"
  print "%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%.4f\t%.4f\t%.4f" % (fa, cor, unc_aligned, cor_aligned, (cor_aligned-unc_aligned), incorrect_uncorrected, correct_uncorrected, incorrect_corrected, correct_corrected, tp, fp, tn, fn, (float(tp)/(tp + fn)), (float(tn)/(tn + fp)), gain)


if __name__ == "__main__":
  parser = argparse.ArgumentParser("Compute true/false positive and gain stats for long-read error corrections")
  parser.add_argument("uncorrected", help="m5 file for original reads")
  parser.add_argument("corrected", help="m5 file for corrected reads")
  parser.add_argument("fasta", help="Uncorrected fasta")
  parser.add_argument("method", help="Method used to correct reads")
  parser.add_argument("--verbose", help="Verbosity (every read)", default=False, action="store_true")
  args = parser.parse_args()
  main(args.uncorrected, args.corrected, args.fasta, args.method, args.verbose)
