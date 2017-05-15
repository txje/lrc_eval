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


Evaluation:
* tef\_stats.py
  * Computes error correction statistics given uncorrected and corrected TEF files
