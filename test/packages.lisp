;; -*- fill-column:102; comment-column:24; -*-

(defpackage :test-cl-sgp4
  (:use :common-lisp :cl-sgp4)
  (:documentation
   "Test Common Lisp implementation of SGP4 satellite propagation algorithm.

Function TEST-SGP4 runs test cases drawn from files in the format used by the programs of
Vallado2006 (see sgp4.lisp for references) and produces a summary table.  The tle and comparison
output files can be specified.  The defaults are the names in Vallado's code, and these default files
are provided in the test directory.

If you use the defaults and hope to find the files in this directory, *default-pathname-defaults*
should point to this directory.  With emacs/slime, for example, this will happen if slime is started
from a buffer containing a file in this directory.

Function TIME-TEST-SGP4 runs the same test cases as TEST-SGP4 but runs the propagation phase n
times (typically 100) and prints the time required for that.

Function TEST-SGP4-1 runs a single case with 2-line elements and times defined by hardcoded values,
perhaps for debugging.

Function README-EXAMPLE does the example in the README file.

Function TLE-FILE-EXAMPLE shows example processing for a complete TLE file.")
  (:export
   #:test-sgp4
   #:time-test-sgp4
   #:test-sgp4-1
   #:readme-example
   #:tle-file-example))
