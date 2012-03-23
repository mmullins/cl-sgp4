;; -*- fill-column:102; comment-column:24; -*-

(in-package :cl-sgp4)


(defstruct tledata
  ;; Data derived from a TLE.  This is raw data in TLE units.
  ;; Note that the identifier line common in tle files is not here.
  ;; It must be managed separately.

  ;; line 1
  (num nil :type fixnum)
  (classification nil :type character)
  (designator-launch-year nil :type fixnum)
  (designator-num-of-year nil :type fixnum)
  (designator-piece-of-launch nil :type string)
  (iyr-norad-epoch nil :type fixnum) ;2-digit year
  (doy-norad-epoch nil :type double-float)
  (rev-per-day2-dmm-dt-over-2 nil :type double-float)   ;first derivative of mean motion
  (rev-per-day3-dmm-dt2-over-6 nil :type double-float)  ;second derivative of mean motion
  (Bstar nil :type double-float)
  (ephermeris-type nil :type character)
  (element-num nil :type fixnum)
  (line1-checksum nil :type fixnum)

  ;; line 2
  ;; num for line2 must be same as for line 1 so is not included separately
  (incl-deg nil :type double-float)
  (ra-of-asc-node-deg nil :type double-float)
  (ecc nil :type double-float)
  (arg-perigee-deg nil :type double-float)
  (mean-anomaly-deg nil :type double-float)
  (mm-rev-per-day nil :type double-float) ;mean motion in rev/day
  (rev-at-epoch nil :type fixnum)
  (line2-checksum nil :type fixnum))

(defun tle-line-calculated-checksum (line)
  "Return the checksum used by TLEs for LINE."
  (mod (loop for i below (min 68 (length line))
          sum (let* ((c (char line i))
                     (digit (digit-char-p c)))
                (cond (digit digit)
                      ((char= c #\-) 1)
                      (t 0))))
       10))

(defun substring-parse-integer (s i-0 i-lim)
  "Parse an integer out of the characters of string S from character index I-0 to just before
character index I-LIM.  Parse blank as 0."
  (let ((s1 (string-trim " " (subseq s i-0 i-lim))))
    (cond ((zerop (length s1)) 0)
          (t (parse-integer s1)))))

(defun substring-parse-double (s i-0 i-lim)
  "Parse a double-float out of the characters of string S from character index I-0 to just before
character index I-LIM."
  (let ((ff *read-default-float-format*))
    (setf *read-default-float-format* 'double-float)
    (let ((d (read-from-string (subseq s i-0 i-lim))))
      (setf *read-default-float-format* ff)
      (the double-float d))))

(defun make-tledata-from-lines (line1 line2)
  "Return a tledata structure containing the data in TLE lines LINE1 and LINE2."
  (when (char/= #\1 (char line1 0))
    (error "card number incorrect on line1: ~s" line1))
  (when (char/= #\2 (char line2 0))
    (error "card number incorrect on line2: ~s" line2))
  (when (not (loop for i across #(1 8 17 32 43 52 61 63) always (char= (char line1 i) #\space)))
    (error "character where space expected on line1: ~s" line1))
  (when (not (loop for i across #(1 7 16 25 33 42 51) always (char= (char line2 i) #\space)))
    (error "character where space expected on line2: ~s" line2))
  (when (not (loop for i across #(23 34) always (char= (char line1 i) #\.)))
    (error "character where '.' expected on line1: ~s" line1))
  (when (not (loop for i across #(11 20 37 46 54) always (char= (char line2 i) #\.)))
    (error "character where '.' expected on line2: ~s" line2))
  ;; (when (not-every? #(= \- (char line1 %)) [50 59])
  ;;   (error "character where '-' expected on line 1"))
  (when (not (loop for i across #(33 44 50 53 59) always
                  (let ((c (char line1 i)))
                    (or (char= c #\-) (char= c #\+) (char= c #\space)))))
    (error "character where '-' '+' or space expected on line1: ~s" line1))
  (let (;; line1
        (num (substring-parse-integer line1 2 7))
        (classification (char line1 7))
        (designator-launch-year (substring-parse-integer line1 9 11))
        (designator-num-of-year (substring-parse-integer line1 11 14))
        (designator-piece-of-launch (string-trim " " (subseq line1 14 17)))
        (iyr-norad-epoch (substring-parse-integer line1 18 20))
        (doy-norad-epoch (substring-parse-double line1 20 32))
        (rev-per-day2-dmm-dt-over-2 (substring-parse-double line1 33 43))
        (rev-per-day3-dmm-dt2-over-6 (* ;format is 12345-x --> 0.12345*10^x
                                      (substring-parse-integer line1 44 50)
                                      (expt 10.0d0 (substring-parse-integer line1 50 52))))
        (Bstar (*       ;format is 12345-x --> 0.12345*10^x
                (* (substring-parse-integer line1 53 59) 1.0d-5)
                (expt 10.0d0 (substring-parse-integer line1 59 61))))
        (ephermeris-type (char line1 62))
        (element-num (substring-parse-integer line1 64 68))
        (line1-checksum (substring-parse-integer line1 68 69))
        ;; line2
        (incl-deg (substring-parse-double line2 8 16))
        (ra-of-asc-node-deg (substring-parse-double line2 17 25))
        (ecc (* 1.0d-7 (coerce (substring-parse-integer line2 26 33) 'double-float)))
        (arg-perigee-deg (substring-parse-double line2 34 42))
        (mean-anomaly-deg (substring-parse-double line2 43 51))
        (mm-rev-per-day (substring-parse-double line2 52 63))
        (rev-at-epoch (substring-parse-integer line2 63 68))
        (line2-checksum (substring-parse-integer line2 68 69)))
    (when (or (< doy-norad-epoch 0.0d0) (>= doy-norad-epoch 366.0d0))
      (error "epoch day-of-year out of range on line1: ~s" line1))
    ;; TODO - other range checks
    (when (/= line1-checksum (tle-line-calculated-checksum line1))
      (error "checksum error for line1: ~s" line1))
    (when (/= line2-checksum (tle-line-calculated-checksum line2))
      (error "checksum error for line2: ~s" line2))
    (make-tledata
     ;; line1
     :num num :classification classification
     :designator-launch-year designator-launch-year
     :designator-num-of-year designator-num-of-year
     :designator-piece-of-launch designator-piece-of-launch
     :iyr-norad-epoch iyr-norad-epoch :doy-norad-epoch doy-norad-epoch
     :rev-per-day2-dmm-dt-over-2 rev-per-day2-dmm-dt-over-2
     :rev-per-day3-dmm-dt2-over-6 rev-per-day3-dmm-dt2-over-6
     :Bstar Bstar
     :ephermeris-type ephermeris-type :element-num element-num :line1-checksum line1-checksum
     ;; line 2
     :incl-deg incl-deg :ra-of-asc-node-deg ra-of-asc-node-deg :ecc ecc
     :arg-perigee-deg arg-perigee-deg :mean-anomaly-deg mean-anomaly-deg
     :mm-rev-per-day mm-rev-per-day :rev-at-epoch rev-at-epoch :line2-checksum line2-checksum)))

(defun tledata-t50-epoch (tledata)
  "Return the epoch time of TLEDATA in the form used by SGP4, i.e., days after 1950.0 ."
  (let* ((iyr-norad-epoch (tledata-iyr-norad-epoch tledata))
         (iyr (if (>= iyr-norad-epoch 50) (+ 1900 iyr-norad-epoch) (+ 2000 iyr-norad-epoch))))
    (+ (ymdhmsZ->t50 iyr 1 1 0 0 0) (1- (tledata-doy-norad-epoch tledata))))
                                 ;;; 1- because doy starts at 1
  )

(defun tledata->t50-v-el0-Bstar (tledata)
  "Return (values t50 v-el0 Bstar) derived from TLEDATA.  Here
  - t50 is the TLE epoch time in days after 1950.0
  - v-el0 is the orbital elements vector at epoch.  These are the so-called Kozai elements, the
    ones directly derived from the TLE
  - Bstar is the TLE B-star drag parameter for the satellite."
  (let* ((iyr-norad-epoch (tledata-iyr-norad-epoch tledata))
         (iyr (if (>= iyr-norad-epoch 50) (+ 1900 iyr-norad-epoch) (+ 2000 iyr-norad-epoch)))
         (t50 (+ (ymdhmsZ->t50 iyr 1 1 0 0 0)
                   (1- (tledata-doy-norad-epoch tledata)))) ;1- because doy starts at 1
         (n0 (/ (* 2.0d0 pi (tledata-mm-rev-per-day tledata)) 1440.0d0))
         (e0 (tledata-ecc tledata))
         (I0 (deg->rdn (tledata-incl-deg tledata)))
         (O0 (deg->rdn (tledata-ra-of-asc-node-deg tledata)))
         (w0 (deg->rdn (tledata-arg-perigee-deg tledata)))
         (M0 (deg->rdn (tledata-mean-anomaly-deg tledata))))
    (values t50 (vector n0 e0 I0 O0 w0 M0) (tledata-Bstar tledata))))

(defun read-tledata-file (path)
  "Read file at PATH which should contain two-line-element data in typical form, i.e., with a
preceding line for each containing an identifier or name for the satellite.  Return an alist of
conses ( <name> . <tledata )."
  (let ((lines (with-open-file (stream path)
                 (loop for line = (read-line stream nil :eof)
                    until (eq line :eof)
                    collect line))))
    (loop for (line0 line1 line2) on lines by #'cdddr
       for id = (string-trim '(#\Space #\Tab #\Linefeed #\Return) line0)
       collect (cons id (make-tledata-from-lines line1 line2)))))
