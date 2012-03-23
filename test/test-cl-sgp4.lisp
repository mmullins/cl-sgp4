;; -*- fill-column:102; comment-column:24; -*-

;; Copyright (c) Mayes Mullins, 2012
;; mmullins@mullinsenterprises.ca

(in-package :test-cl-sgp4)

(declaim (optimize (debug 3)))


(defun read-tledata.times (path)
  "Read file at PATH, which should be in the form used by the programs of Vallado2006, and return a
list of conses, each a pair of tledata and a list of times at which state is desired for each
satellite."
  (let ((lines (with-open-file (stream path)
                 (loop for line = (read-line stream nil :eof)
                    until (eq line :eof)
                    collect line)))
        (tledata.times '())
        (ff *read-default-float-format*))
    (setf *read-default-float-format* 'double-float)
    (loop until (endp lines) do
         (let ((line1 (first lines)))
           (if (char= (char line1 0) #\#)
               (setf lines (cdr lines)) ;comment line
               (let* ((line2+times (second lines))
                      (line2 (subseq line2+times 0 69))
                      (tledata (make-tledata-from-lines line1 line2))
                      (time-string (subseq line2+times 69))
                      (times
                       (multiple-value-bind (from-time next-i-1)
                           (read-from-string time-string t nil :start 0)
                         (multiple-value-bind (to-time next-i-2)
                             (read-from-string time-string t nil :start next-i-1)
                           (multiple-value-bind (by-time)
                               (read-from-string time-string t nil :start next-i-2)
                             (let* ((times-1 (loop for time from from-time to to-time by by-time
                                                collect (coerce time 'double-float)))
                                    ;; make sure first time is 0 and last is to-time since comparison
                                    ;; program does that
                                    (times-2 (if (= (car (last times-1)) to-time)
                                                 times-1
                                                 (append times-1 (list to-time))))
                                    (times-3 (if (zerop (first times-2))
                                                 times-2
                                                 (cons 0.0d0 times-2))))
                               times-3))))))
                 (push (cons tledata times) tledata.times)
                 (setf lines (cddr lines))))))
    (setf *read-default-float-format* ff)
    (nreverse tledata.times)))

(defun read-doubles-from-string (s)
  "Return a list of doubles read from string S."
  (let ((ff *read-default-float-format*))
    (setf *read-default-float-format* 'double-float)
    (let* ((forms       ;list of top-level forms in the file
            (let ((eof (gensym)))
              (with-input-from-string (s-stream s)
                (loop for form = (read-preserving-whitespace s-stream nil eof)
                        ;probably just read is enough here
                   until (eq form eof)
                   collect form))))
           (doubles (mapcar (lambda (form) (coerce form 'double-float)) forms)))
      (setf *read-default-float-format* ff)
      doubles)))

(defun read-sgp4-outputs (path)
  "Read a file at PATH, which should be in the form used by the programs of Vallado2006,
i.e. is in the format:
{
 one line with a satellite number in the first 20 characters and nothing after
 { lines each with time, 3 positions, 3 velocities
   (the total line taking more than 20 characters)
   }*
 }*
Return a list of conses, each a cons/pair of a num and a cons, with the second
cons ( <time> . <6-element-state-vector> )."
  (let ((lines (with-open-file (stream path)
                 (loop for line = (read-line stream nil :eof)
                    until (eq line :eof)
                    collect line))))
    (loop for (line1 . rest) = lines
       while line1
       for num = (read-from-string line1)
       for l-cons = (loop for line = (car rest)
                       while (> (length line) 20)
                       for l7 = (read-doubles-from-string line)
                       for time = (car l7)
                       for v-state = (map 'vector #'identity (cdr l7))
                       collect (cons time v-state)
                       do (setf rest (cdr rest)))
       for v-cons = (map 'vector #'identity l-cons)
       collect (cons num v-cons)
       do
         (setf lines rest)
         (when (> (length line1) 20)
           (error "expected satellite number line but this seems too long: ~s~%" line1))
         )))

(defun test-sgp4 (&key
                  (spec-path #P"SGP4-VER.TLE")
                  (compare-path #P"tforver.out")
                  (delta-limit 1.0d-3) (delta-v-limit 1.0d-6))
  "Run test cases in the file at path SPEC-PATH and compare the results with those in the file at path
COMPARE-PATH.  These are in the form used and produced by the test programs of Vallado2006.  Mark with
asterisks those values of delta greater than DELTA-LIMIT meters and of delta-v greater than
DELTA-V-LIMIT meters/second."
  (let* ((tledata.times (read-tledata.times spec-path))
         (sgp4datas (loop for (tledata . times) in tledata.times ;takes 0.001 sec
                       ;; (format t "~a~%" (tledata-num tledata))
                       collect
                         (cons (tledata-num tledata)
                               (multiple-value-bind (t50 v-el0 Bstar)
                                   (tledata->t50-v-el0-Bstar tledata)
                                 (let* ((v-grav (make-v-grav-wgs-72))
                                        (sgp4data (init-sgp4 v-grav t50 v-el0 Bstar)))
                                   sgp4data)))))
         (outputs (time (loop for (tledata . times) in tledata.times
                           for (nil . sgp4data) in sgp4datas
                           ;; (format t "~a~%" (tledata-num tledata))
                           collect
                             (cons
                              (tledata-num tledata)
                              (map 'vector #'identity
                                   (loop for time in times
                                      collect
                                        (cons
                                         time
                                         (when sgp4data
                                           (handler-case
                                               (->v-state-km-s
                                                (propagate-sgp4 sgp4data time))
                                             (altitude-lt-0-error ()
                                               (format t "altitude < 0 for num = ~a at time = ~a~%"
                                                       (tledata-num tledata) time)
                                               nil))))))))))
         (comparison-outputs (read-sgp4-outputs compare-path)))
    (format t "~10@a  ~12@a ~1@a ~1@a ~1@a ~1@a ~6@a / ~6@a / ~6@a ~14@a ~1a ~14@a ~1a~%"
            "num" "period [min]" "D" "R" "T" "L" "# tried" "# done" "# comp"
            "delta [m]" " " "  delta-v [m/s]" " ")
    (loop
       for max-delta = 0.0d0
       for max-delta-v = 0.0d0
       for (tledata . nil) in tledata.times ;used for output labelling
       for (nil . sgp4data) in sgp4datas    ;used for output labelling
       for (num . v-conses) in outputs
       for (comparison-num . comparison-v-conses) in comparison-outputs do
         (when (/= num comparison-num) (error "nums ~a and ~a don't match~%" num comparison-num))
         (loop
            for (time . nil) across v-conses
            for (comparison-time . nil) across comparison-v-conses
            do (when (> (abs (- time comparison-time)) 1.0d-6)
                 (error "times ~a ~a not equal~%" time comparison-time)))
         (loop
            for (time . v-state) across v-conses
            while v-state
            for (comparison-time . comparison-v-state) across comparison-v-conses
            for x = (aref v-state 0)
            for y = (aref v-state 1)
            for z = (aref v-state 2)
            for comparison-x = (aref comparison-v-state 0)
            for comparison-y = (aref comparison-v-state 1)
            for comparison-z = (aref comparison-v-state 2) do
              (let ((max-1 (max (abs (- x comparison-x))
                                (abs (- y comparison-y))
                                (abs (- z comparison-z)))))
                (when (> max-1 max-delta)
                  (setf max-delta max-1))))
         (loop
            for (time . v-state) across v-conses
            while v-state
            for (comparison-time . comparison-v-state) across comparison-v-conses
            for vx = (aref v-state 3)
            for vy = (aref v-state 4)
            for vz = (aref v-state 5)
            for comparison-vx = (aref comparison-v-state 3)
            for comparison-vy = (aref comparison-v-state 4)
            for comparison-vz = (aref comparison-v-state 5) do
              (let ((max-1 (max (abs (- vx comparison-vx))
                                (abs (- vy comparison-vy))
                                (abs (- vz comparison-vz)))))
                (when (> max-1 max-delta-v)
                  (setf max-delta-v max-1))))
         (let ((period-min (/ 1440.0 (tledata-mm-rev-per-day tledata)))
               (deep-space? (and sgp4data (sgp4data-deep-space? sgp4data)))
               (resonance (and sgp4data (sgp4data-resonance sgp4data)))
               (drop-terms? (and sgp4data (sgp4data-drop-terms? sgp4data)))
               (near-lyddane? (< (* pi (/ (tledata-incl-deg tledata) 180.0d0)) 0.22d0))
               (v-conses-done (cond ((null sgp4data) 0)
                                    (t (length (remove (lambda (c) (null (cdr c))) v-conses))))))
           (format t "~10d  ~12,4f ~a ~a ~a ~a ~6d / ~6d / ~6d ~14,2e ~a ~14,2e ~a~%"
                   num period-min
                   (if deep-space? "D" " ")
                   (case resonance (:synchronous "S") (:half-day "H") (otherwise " "))
                   (if drop-terms? "T" " ")
                   (if near-lyddane? "L" " ")
                   (length v-conses) v-conses-done (length comparison-v-conses)
                   (* 1000.0d0 max-delta)
                   (if (> (* 1000.0d0 max-delta) delta-limit) "*" " ")
                   (* 1000.0d0 max-delta-v)
                   (if (> (* 1000.0d0 max-delta-v) delta-v-limit) "*" " "))))))

(defun time-test-sgp4 (&key
                       (spec-path #P"SGP4-VER.TLE")
                       (n 100))
  (declare (optimize (speed 3) (safety 0)))
  "Apply cl:time to N repititions of the test cases in the file at path SPEC-PATH."
  (let* ((tledata.times (read-tledata.times spec-path))
         (sgp4datas (loop for (tledata . times) in tledata.times ;takes 0.001 sec
                       ;; (format t "~a~%" (tledata-num tledata))
                       collect
                         (cons (tledata-num tledata)
                               (multiple-value-bind (t50 v-el0 Bstar)
                                   (tledata->t50-v-el0-Bstar tledata)
                                 (let* ((v-grav (make-v-grav-wgs-72))
                                        (sgp4data (init-sgp4 v-grav t50 v-el0 Bstar)))
                                   sgp4data))))))
    (format t "doing ~a calls to propagate~%"
            (* n (loop for (nil . times) in tledata.times sum (length times))))
    (time
     (loop for (tledata . times) in tledata.times
        for (nil . sgp4data) in sgp4datas
        collect
        (cons
         (tledata-num tledata)
         (map 'vector #'identity
              (loop for time in times
                 collect
                 (cons
                  time
                  (when sgp4data
                    (handler-case
                        (->v-state-km-s
                         (propagate-sgp4 sgp4data time))
                      (altitude-lt-0-error ()
                        nil))
                    (dotimes (i n)
                      (handler-case
                          (->v-state-km-s
                           (propagate-sgp4 sgp4data time))
                        (altitude-lt-0-error ()
                          nil)))))))))))
  (values))

(defun test-sgp4-1 (&key
                    (two-lines
                     ;; '("1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753"
                     ;;   "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667")

                     ;; '("1 22312U 93002D   06094.46235912  .99999999  81888-5  49949-3 0  3953"
                     ;;   "2 22312  62.1486  77.4698 0308723 267.9229  88.7392 15.95744531 98783")

                     '("1 23333U 94071A   94305.49999999 -.00172956  26967-3  10000-3 0    15"
                       "2 23333  28.7490   2.3720 9728298  30.4360   1.3500  0.07309491    70")

                     ;; '("1 23599U 95029B   06171.76535463  .00085586  12891-6  12956-2 0  2905"
                     ;;   "2 23599   6.9327   0.2849 5782022 274.4436  25.2425  4.47796565123555")

                     ;; '("1 28350U 04020A   06167.21788666  .16154492  76267-5  18678-3 0  8894"
                     ;;   "2 28350  64.9977 345.6130 0024870 260.7578  99.9590 16.47856722116490")

                     ;; '("1 28872U 05037B   05333.02012661  .25992681  00000-0  24476-3 0  1534"
                     ;;   "2 28872  96.4736 157.9986 0303955 244.0492 110.6523 16.46015938 10708")

                     ;; '("1 29141U 85108AA  06170.26783845  .99999999  00000-0  13519-0 0   718"
                     ;;   "2 29141  82.4288 273.4882 0015848 277.2124  83.9133 15.93343074  6828")
                     )
                    (from-time 0.0d0) (to-time 0.0d0) (by-time 5.0d0))
  "Run a single case with 2-line elements defined by hard-coded line1 and line2."
  (let* ((line1 (first two-lines))
         (line2 (second two-lines))
         (tledata (make-tledata-from-lines line1 line2)))
    (format t "line1 = ~s~%line2 = ~s~%tledata =~%~a~%" line1 line2 tledata)
    (multiple-value-bind (t50 v-el0 Bstar) (tledata->t50-v-el0-Bstar tledata)
      (format t "t50 = ~a~%v-el0 = ~a~%" t50 v-el0)
      (let* ((v-grav (make-v-grav-wgs-72))
             (sgp4data (init-sgp4 v-grav t50 v-el0 Bstar)))
        (format t "sgp4data = ~%~a~%" sgp4data)
        (loop
           for time from from-time to to-time by by-time
           for v-state = (propagate-sgp4 sgp4data time)
           for v-location-km = (->v-location-km v-state)
           for v-velocity-km/s = (->v-velocity-km/s v-state) do
           (format t "~16,8f" time)
           (loop for i from 0 below 3 do
                (format t "~17,8f" (aref v-location-km i)))
           (loop for i from 0 below 3 do
                (format t "~15,9f" (aref v-velocity-km/s i)))
           (format t "~%"))))))

(defun readme-example ()
  "Run the example in the README file."
  (let ((tledata
         (make-tledata-from-lines
          "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753"
          "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667"))
        (v-grav (make-v-grav-wgs-72)))
    (multiple-value-bind (t50 v-el0 Bstar) (tledata->t50-v-el0-Bstar tledata)
      (let ((sgp4data (init-sgp4 v-grav t50 v-el0 Bstar)))
        (format t "~10@a  ~10@a~10@a~10@a  ~10@a~10@a~10@a~%"
                "TIME [min]" "X [km]" "Y [km]" "Z [km]"
                "VX [km/s]" "VY [km/s]" "VZ [km/s]")
        (loop
           for time from 0.0d0 to 360.0d0 by 60.0d0
           for v-state-km-s = (->v-state-km-s (propagate-sgp4 sgp4data time))
           do (format t "~10,3f  " time)
             (loop for i from 0 below 3
                do (format t "~10,3f" (aref v-state-km-s i)))
             (format t "  ")
             (loop for i from 3 below 6
                do (format t "~10,6f" (aref v-state-km-s i)))
             (format t "~%"))))))

(defun tle-file-example (&key
                         (path #P"visual.txt")
                         ;; can download from http://www.celestrak.com/NORAD/elements/visual.txt
                         (tle-alist (read-tledata-file path))
                         (ids '(25544))
                         (from-universal-time (get-universal-time))
                         (to-universal-time (+ from-universal-time 3600.0d0))
                         (by-seconds 300.0d0))
  "Take TLE data from either the alist TLE-ALIST, which should map id to tledata structures,
or read it from the TLE data file at PATH.  In the latter case, the file should be in the usual format
in which TLE data is provided, i.e., with 3 lines per object, the 'zero-th' line containing an
identifying name for the object, and the other 2 lines containing the TLE.  If both are provided, use
only TLE-ALIST and ignore PATH.

Out of this TLE data, process the objects identified by the list IDS.  Each object is identified by
either an integer which is the identifying number of the object in its TLE, or a string which is the
key in TLE-ALIST for the object.  As a special case, if IDS is just a single integer or string,
process just that one object.

For each object to be processed, get its state at times from FROM-UNIVERSAL-TIME, to
TO-UNIVERSAL-TIME, with steps BY-SECONDS.  The universal times are those returned by such Common Lisp
functions as CL:GET-UNIVERSAL-TIME and CL:ENCODE-UNIVERSAL-TIME.  For each object and each time, print
- the date and time
- the time minus the epoch time of the TLE in minutes
- the object's location in ECI coordinates
- the object's location as longitude, positive east; latitude, positive north, and altitude."
  (let* ((ids (if (listp ids) ids (list ids)))
         (tledatas (map 'list (lambda (id)
                                (or (and (stringp id) (cdr (assoc id tle-alist :test #'string=)))
                                    (and (integerp id) (cdr (rassoc id tle-alist :key #'tledata-num)))))
                        ids))
         (v-grav (make-v-grav-wgs-72)))
    (loop for id in ids
       for tledata in tledatas
       for sgp4data = (multiple-value-bind (t50 v-el0 Bstar) (tledata->t50-v-el0-Bstar tledata)
                        (init-sgp4 v-grav t50 v-el0 Bstar))
       for epoch-universal-time = (t50->universal-time (sgp4data-t50 sgp4data))
       do
       (format t "id = ~a~%" id)
       (format t "~19@a  ~12@a  ~12@a~12@a~12@a  ~12@a~12@a~12@a~%"
               "DATE/TIME" "-EPOCH [min]" "X [km]" "Y [km]" "Z [km]"
               "ELONG [DEG]" "NLAT [DEG]" "ALT [KM]")
       (loop
          for universal-time from from-universal-time to to-universal-time by by-seconds
          for minutes-after-epoch = (/ (- universal-time epoch-universal-time) 60.0d0)
          for v-state = (propagate-sgp4 sgp4data minutes-after-epoch)
          for v-location-km = (->v-location-km v-state)
          do
          (multiple-value-bind (sec min hr date month year)
              (decode-universal-time (floor universal-time))
            (format t "~4d-~2,'0d-~2,'0d ~2,'0d:~2,'0d:~2,'0d  " year month date hr min sec))
          (format t "~12,3f  " minutes-after-epoch)
          (loop for i from 0 below 3 do
               (format t "~12,3f" (aref v-location-km i)))
          (format t "  ")
          (multiple-value-bind (elong lat alt)
              (txyz->lla universal-time
                         (aref v-location-km 0) (aref v-location-km 1) (aref v-location-km 2))
            (format t "~12,3f~12,3f~12,3f" (rdn->deg elong) (rdn->deg lat) alt))
          (format t "~%")))))
