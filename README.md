A few tools to evaluate error correction of long reads (Pacbio, nanopore), largely mirroring the [Error Correction Evaluation Toolkit](http://aluru-sun.ece.iastate.edu/doku.php?id=ecr)


From the ECET paper:

> We use the following measures for each program:
> number of erroneous bases identified and
> successfully corrected (true positives, TP), correct
> bases wrongly identified as errors and changed
> (false positives, FP), and erroneous bases that were
> either uncorrected or falsely corrected (false negatives,
> FN). We report sensitivity and specificity for
> each program. Then, we combine these into the gain
> metric [21], defined by gain = (TP - FP) /
> (TP + FN), which is the percentage of errors
> removed from the data set by the error-correction
> program. A negative gain value indicates that more
> errors have been introduced due to false corrections,
> which is not captured by measures such as sensitivity
> and specificity.


Utilities:
* maf2tef.py
  * converts MAF to TEF format
* sam2tef.py
  * converts SAM to TEF format
* m52tef.py
  * converts BLASR -m5 format to TEF format
* remap_m5.py
  * rewrites the read names in a FASTA file according to the renaming scheme for several long read error correction methods
  * it's easier to compare post- to pre-corrected sequences if the names are consistent...

Plumbing:
* fasta.py
  * A very simple FASTA file API
* aln_formats.py
  * Provides a common API to parse and iterate through alignment formats, including MAF, m4, and m5

Statistics can be computed directly from several alignment formats, with slightly different capabilities:
* tef_stats.py
  * Computes error correction statistics given uncorrected and corrected TEF files, in line with original ECET
* maf_stats.py
* m5_stats.py
  * THIS IS THE RECOMMENDED METHOD and the method used in the FMLRC paper
  * Statistics are computed directly from -m5 format, allowing BLASR results to be used directly and loci compared relative to the reference sequence
