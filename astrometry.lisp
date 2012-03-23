;; -*- fill-column:102; comment-column:24; -*-

(in-package :cl-sgp4)

;; (declaim (optimize (debug 3)))
(declaim (optimize (speed 3) (safety 0)))


;;;; ANGLE CONVERSIONS

(defun deg->rdn (deg)
  "Convert DEG degrees to radians."
  (* deg (/ pi 180.0d0)))

(defun rdn->deg (rdn)
  "Convert RDN radians to degrees."
  (* rdn (/ 180.0d0 pi)))

(declaim (inline deg->rdn rdn->deg))

(defun dms->rdn (deg arcmin arcsec)
  "Convert degrees, arc-minutes and arc-seconds to radians."
  (deg->rdn (+ deg (/ arcmin 60.0) (/ arcsec 3600.0))))

(defun ndms->rdn (sign deg arcmin arcsec)
  "Convert degrees, arc-minutes and arc-seconds to radians.  The output is given the sign of SIGN."
  (let ((rdn (dms->rdn deg arcmin arcsec)))
    (if (< sign 0) (- rdn) rdn)))

(defun rdn->dms (rdn)
  "Return 3 values:
  - degrees as an integer
  - arc-minutes as an integer, 0..59
  - arc-seconds as a double, [0.0 .. 60.0)
corresponding to an angle given by the absolute value of RDN radians."
  (let* ((deg (rdn->deg (abs rdn)))
         (whole-deg (floor deg))
         (min (* 60.0 (- deg whole-deg)))
         (whole-min (floor min))
         (sec (* 60.0 (- min whole-min))))
    (values whole-deg whole-min sec)))

(defun rdn->ndms (rdn)
  "Return a 4 values:
  - sign as -1 or +1
  - degrees as an integer
  - arc-minutes as an integer, 0..59
  - arc-seconds as a double, [0.0 .. 60.0)
corresponding to an angle RDN radians."
  (let* ((sgn (if (< rdn 0) -1 +1))
         (deg (rdn->deg (abs rdn)))
         (whole-deg (floor deg))
         (min (* 60.0 (- deg whole-deg)))
         (whole-min (floor min))
         (sec (* 60.0 (- min whole-min))))
    (values sgn whole-deg whole-min sec)))

(defun rdn-0-2pi (rdn)
  "Reduce angle RDN radians to the range [0..2*pi), i.e., greater than or equal to zero, less than
2*pi.  Return a double-float."
  (let ((r (rem rdn (* 2.0d0 pi))))
    (if (>= r 0.0) r (+ r (* 2.0d0 pi)))))

(defun rdn-mpi-ppi (rdn)
  "Reduce angle RDN radians to the range (-pi..+pi], i.e., greater than -pi, less than or equal to
+pi.  Return a double-float."
  (let ((r (rem rdn (* 2.0d0 pi))))
    (if (<= (abs r) pi)
        (if (> r (- pi)) r pi)
        (if (< rdn 0.0d0) (+ r (* 2.0d0 pi)) (- r (* 2.0d0 pi))))))

(declaim (inline rdn-0-2pi rdn-mpi-ppi))


;;;; SPHERICAL GEOMETRY

(defun v->phi-theta (v)
  "Return (values phi theta), spherical coordinates corresponding to the direction of 3-vector V.
Here phi is longitude or right ascension, theta is latitude or declination (not co-latitude as in
usual spherical coordinates)."
  (let* ((x (aref v 0))
         (y (aref v 1))
         (z (aref v 2))
         (r-xy (sqrt (+ (* x x) (* y y)))))
    (if (zerop r-xy)
        (values 0.0 (* 0.5d0 pi)) ;phi is arbitrary at this singularity
        (values (atan y x) (atan (/ z r-xy))))))

(defun normalized-v (v)
  "Return a vector in the same direction as input vector V, but normalized, i.e., a unit
vector in the same direction as V."
  (let* ((n (length v))
         (mag (sqrt (loop for i below n
                       sum (expt (aref v i) 2))))
         (v1 (make-array n :element-type (array-element-type v))))
    (dotimes (i n) (setf (aref v1 i) (/ (aref v i) mag)))
    v1))

(defun phi-theta->normalized-v (phi theta)
  "Return a 3-vector containing a normalized vector in the direction of spherical coordinates PHI,
THETA.  Here PHI is longitude or right ascension, THETA is latitude or declination (not co-latitude as
in usual spherical coordinates)."
  (let ((cos-theta (cos theta)))
    (vector (* cos-theta (cos phi)) (* cos-theta (sin phi)) (sin theta))))

(defun rdn-between-normalized-v (v-nrm1 v-nrm2)
  "Return a double-float giving the angle in radians between the two normalized 3-vectors V-NRM1 and
V-NRM2.  Normalization is not checked, and the value returned will be in error if the vectors are not
normalized.  If in doubt, use rdn-between-v."
  (let ((mod-sq-half (* 0.25 (loop for i below 3 sum
                                  (expt (- (aref v-nrm1 i) (aref v-nrm2 i)) 2)))))
    (atan (sqrt mod-sq-half)
          (sqrt (max 0.0d0 (- 1.0d0 mod-sq-half))))))

(defun rdn-between-v (v-1 v-2)
  "Return a double-float giving the angle in radians between the two 3-vectors V-1 and V-2."
    (rdn-between-normalized-v (normalized-v v-1) (normalized-v v-2)))

(defun rdn-between-phi-theta (phi1 theta1 phi2 theta2)
  "Return a double-float giving the angle in radians between the two directions specified by the
spherical coordinate pairs (PHI1, THETA1) and (PHI2,THETA2).  Here the PHI-s are longitude or right
ascension, the THETA-s latitude or declination (not co-latitude as in usual spherical coordinates)."
    (rdn-between-normalized-v (phi-theta->normalized-v phi1 theta1)
                              (phi-theta->normalized-v phi2 theta2)))


;;;; TIME CONVERSIONS

;; iyr = integer 4-digit year
;; imoy = integer month of year, 1..12
;; idom = integer day of month, 1..28,29,30,31 as appropriate for the month and year
;; ihod = integer hour of day, 0..23
;; imoh = integer minute of hour, 0..59
;; som  = second of minute, [0.0 .. 60.0)

;; Note that the order of these in parameter lists and return values is reversed from that of the
;; Common Lisp functions CL:DECODE-UNIVERSAL-TIME, CL:GET-DECODED-TIME and CL:ENCODE-UNIVERSAL-TIME,
;; but is probably more consistent with conventional use.

;; JD = Julian date

(defun leap-year-p (iyr)
  "Return T if IYR is a leap year, nil if not."
  (let ((iyr (floor iyr)))
    (and (zerop (mod iyr 4))
         (or (not (zerop (mod iyr 100)))
             (zerop (mod iyr 400))))))

(let ((nday-in-imoy-not-ly #(0 31 28 31 30 31 30 31 31 30 31 30 31)))
                                        ;note - offset by 1
  (defun nday-in-iyr-imoy (iyr imoy)
    "Return the number of days in month IMOY of year IYR."
    (when (or (< imoy 1) (> imoy 12)) (error "imoy ~a out of range, should be 1..12" imoy))
    (cond ((/= imoy 2) (aref nday-in-imoy-not-ly imoy))
          ((leap-year-p iyr) 29)
          (t 28))))

(defun ymd->iJD (iyr imoy idom)
  "Return integer Julian day starting at Greenwich noon of Gregorian calendar date IYR (year,
1582...), IMOY (month number, 1..12) IDOM (day of month, 1..28,29,30,31 as appropriate for the month).
Actually any input that gives a Julian day greater than or equal to 0 is acceptable with JD = 0
corresponding to Greenwich noon of 1 Jan 4713 BPE.  However, the Gregorian calendar officially began
on Fri 15 Oct 1582, 1 day after Thurs 4 Oct 1528 Julian calendar.  Switch took place in Germany in
1698, England in 1752, for example."
  ;; T C Van Flandern, K F Pulkkinen
  ;; "Low-Precision Formulae for Planetary Positions"
  ;; Astrophys J Supplement Series v41, p391-411, Nov79
  ;; epoch time is in days after Jan 0, 1950 (1949-12-31T00:00:00.000Z) = Julian day 2433281.5
  (when (or (< imoy 1) (> imoy 12)) (error "imoy ~a out of range, should be 1..12" imoy))
  (when (or (< idom 1) (> idom (nday-in-iyr-imoy iyr imoy)))
    (error "idom ~a out of range, should be 1..~a" idom (nday-in-iyr-imoy iyr imoy)))
  (let (;; make sure inputs are integers, which is necessary for the following algorithm
        (iyr (floor iyr))
        (imoy (floor imoy))
        (idom (floor idom)))
    (- (+ (* 367 iyr) (floor (* 275 imoy) 9) idom 1721029)
       (floor (* 7 (+ iyr (floor (+ imoy 9) 12))) 4)
       (floor (* 3 (+ 1 (floor (+ iyr (floor (- imoy 9) 7)) 100))) 4))))

(defun iJD->ymd (iJD)
  "Return (values iyr imoy idom) containing integers year, month, day-of-month, corresponding to
Julian day 'IJD'-0.5 to 'IJD'+0.5."
  ;; Henry F Fliegel, Thomas C Van Flandern
  ;; "A Machine Algorithm for Processing Calendar Dates"
  ;; Comm ACM v11, n10, p657, Oct68
  (let* ((iJD (floor iJD))
         (l (+ iJD 68569))
         (n (floor (* 4 l) 146097))
         (l (- l (floor (+ 3 (* 146097 n)) 4)))
         (i (floor (* 4000 (+ l 1)) 1461001))
         (l (+ l (- (floor (* 1461 i) 4)) 31))
         (j (floor (* 80 l) 2447))
         (k (- l (floor (* 2447 j) 80)))
         (l (floor j 11))
         (j (+ j 2 (* -12 l)))
         (i (+ (* 100 (- n 49)) i l)))
    (values i j k)))

(defun hms->day (hr min sec)
  "Return a double representing the number of days corresponding to HR, MIN, SEC, hours, minutes, and
seconds.  Normally MIN is minute-of-hour and SEC is second-of-minute, so both are 0..59.  However,
this isn't essential; the returned value corresponds to the sum of the total times of HR, MIN and SEC."
  (/ (+ (* 3600.0d0 hr) (* 60.0d0 min) sec) 86400.0d0))

(defun day->hms (day)
  "Return 3 values:
  - hours as an integer
  - minutes as an integer, 0..59
  - seconds as a double-float, [0.0 .. 60.0)
corresponding to a duration absolute value of DAY days."
  (let* ((hr (* 24.0 (abs day)))
         (whole-hr (floor hr))
         (min (* 60.0 (- hr whole-hr)))
         (whole-min (floor min))
         (som (* 60.0 (- min whole-min))))
    (values whole-hr whole-min som)))

(defun day->nhms (day)
  "Return 4 values:
  - sign as -1 or +1
  - hours as an integer
  - minutes as an integer, 0..59
  - seconds as a double-float, [0.0 .. 60.0)
corresponding to a duration DAY days."
  (let* ((sign (if (minusp day) -1 +1))
         (hr (* 24.0 (abs day)))
         (whole-hr (floor hr))
         (min (* 60.0 (- hr whole-hr)))
         (whole-min (floor min))
         (som (* 60.0 (- min whole-min))))
    (values sign whole-hr whole-min som)))

(defun ymdhmsZ->JD (iyr imoy idom ihod imoh som)
  "Return a double-float Julian day corresponding to
IYR, IMOY, IDOM (year month-of-year, day-of-month),
IHOD, IMOH, SOM (hour-of-day minute-of-hour, second-of-minute) UT.
All inputs should be integers except for SOM which can be a real."
  (let ((iJD (ymd->iJD iyr imoy idom))
        (day-frac (hms->day ihod imoh som)))
    (+ iJD day-frac -0.5d0)
    ;; - 0.5 adjusts for offset between JD start of 12:00 UT and
    ;; normal day start of 00:00 UT
    ))

(defun JD->ymdhmsZ (JD)
  "Return (values iyr imoy idom ihod imoh som), corresponding to the GMT of Julian date JD.
JD can be an integer or real."
  (when (< JD 0.0d0) (error "JD = ~a but must be >= 0" JD))
  (let* ((JDZ (+ JD 0.5))
         (whole-day (floor JDZ))
         (frac-day (- JDZ whole-day)))
    (multiple-value-bind (iyr imoy idom) (iJD->ymd whole-day)
      (multiple-value-bind (ihod imoh som) (day->hms frac-day)
        (values iyr imoy idom ihod imoh som)))))

;; SGP4 time

(defun t50->JD (t50)
  "Convert SGP4 date value T50 (days after 1950.0) to a Julian date."
  (+ t50 2433281.5d0))

(defun JD->t50 (JD)
  "Convert Julian date JD to an SGP4 date value (days after 1950.0)."
  (- JD 2433281.5d0))

(declaim (inline t50->JD JD->t50))

(defun ymdhmsZ->t50 (iyr imoy idom ihod imoh som)
  "Return a double-float SGP4 date (days after 1950.0) corresponding to
IYR, IMOY, IDOM (year month-of-year, day-of-month),
IHOD, IMOH, SOM (hour-of-day minute-of-hour, second-of-minute) UT.
All inputs should be integers except for SOM which can be a real."
  ;; epoch time is in days after Jan 0, 1950 = Julian day 2433281.5
  ;;     JD 2433281.5 --> 1949-12-31T00:00:00.000Z
  (- (ymdhmsZ->JD iyr imoy idom ihod imoh som) 2433281.5))

(defun t50->ymdhmsZ (t50)
  "Return (values iyr imoy idom ihod imoh som), corresponding to SGP4 date T50 (days after 1950.0)."
  (JD->ymdhmsZ (t50->JD t50)))

;; Common Lisp time

(defun universal-time->JD (u)
  "Return Julian date as a double-float corresponding to Common Lisp 'universal time' U.
From the Common Lisp hyperspec, 'universal time' is 'the number of seconds since midnight, January 1,
1900 GMT (ignoring leap seconds).  As used by normal Common Lisp functions, universal time is a
positive integer.  Here, however, it can be a real for sub-second precision, and can be negative for
pre-1900 times.  U can be produced by CL:GET-UNIVERSAL-TIME and CL:ENCODE-UNIVERSAL-TIME.
Note, however, that the order of parameters in the latter Common Lisp function
 (year, month, day, etc) is reversed from that used by the functions in this package."
  (+ 2415020.5d0 (/ u 86400.0d0)))

(defun JD->universal-time (JD)
  "Return a Common Lisp 'universal time' corresponding to integer or real Julian data JD.
The return value is a double-float.  It can be used by CL:DECODE-UNIVERSAL-TIME if first run through
CL:FLOOR.  Note, however, that the order of return values from that Common Lisp function (year, month,
day, etc) is reversed from that used by the functions in this package."
  (* 86400.0d0 (- JD 2415020.5d0)))

(defun universal-time->t50 (u)
  (JD->t50 (universal-time->JD u)))

(defun t50->universal-time (t50)
  (JD->universal-time (t50->JD t50)))

;; For the equivalent of ymdhmsZ and local time ymdhms conversions to/from universal-time, use
;; Common Lisp functions CL:ENCODE-UNIVERSAL-TIME and CL:DECODE-UNIVERSAL-TIME.
;; NOTE, however, that the order of arguments/values for these is reversed from that used here.


;;;; SIDEREAL TIME

;; Note that the functions below ignore subtleties like precession of the equinoxes and polar wander.

(defun phi-IAU-at-JD (JD)
  "Return GMST (Greenwich Mean Sidereal Time) or GHA (Greenwich Hour Angle), [0..2*pi) using IAU
calculation at time given by JD, Julian date UT1.  Note that this ignores subtleties like precession
of the equinoxes and polar wander."
  (let* ((tut1 (/ (- JD 2451545.0d0) 36525.0d0))
                        ;Julian centuries after jan 1, 2000 12 h epoch (UT1)
         (phi-sec (+ 67310.54841d0 ;GMST or Greenwich Hour Angle (GHA) in seconds of time
                     (* tut1 3.164400184812866d9)
                        ;3.16...d9 = (+ 8640184.812866 (* 876600.0 3600.0))
                     (* tut1 tut1 0.093104d0)
                     (* tut1 tut1 tut1 -6.2d-6))))
    ;; 1 second of time = 360/86400 = 1/240 degree of rotation = 1/43200 of pi radians
    (rdn-0-2pi (* phi-sec (/ pi 43200.0d0)))))

(defun phi-sgp4-at-t50 (t50)
  "Return GMST (Greenwich Mean Sidereal Time) or GHA (Greenwich Hour Angle), [0..2*pi) using SGP4
calculation (1970-Jan-00 reference) at time given by T50 days after Jan 0, 1950. 0 hr. the SGP4
epoch."
  ;; epoch time is in days after Jan 0, 1950 = Julian date 2433281.5
  ;;     (jt<-JD 2433281.5) --> #<DateTime 1949-12-31T00:00:00.000Z>
  (let* ((ts70 (- t50 7305.0d0))
         (ids70 (floor (+ ts70 1.0d-8)))
                        ;1.0d-8 offset is in the canonical SGP4 Fortran program but not in Vallado2006
                        ;(see notes in sgp4.lisp)
         (tfrac (- ts70 ids70))
         (c1 1.72027916940703639d-2)
         (thgr70 1.7321343856509374d0)
         (fk5r 5.07551419432269442d-15)
         (c1p2p (+ c1 (* 2.0d0 pi))))
    (rdn-0-2pi (+ thgr70 (* c1 ids70) (* c1p2p tfrac) (* ts70 ts70 fk5r)))))

(defun phi-IAU-at-t50 (t50)
  "Return GMST (Greenwich Mean Sidereal Time) or GHA (Greenwich Hour Angle), [0..2*pi) using IAU
calculation at time given by T50 days after Jan 0, 1950. 0 hr.  Note that this ignores subtleties
like prcession of the equinoxes and polar wander."
  ;; epoch time t50 is in days after Jan 0, 1950 = Julian date 2433281.5
  (phi-IAU-at-JD (t50->JD t50)))

(defun phi-IAU-at-universal-time (u)
  "Return GMST (Greenwich Mean Sidereal Time) or GHA (Greenwich Hour Angle), [0..2*pi) using IAU
calculation at Common Lisp 'universal time' U."
  (phi-IAU-at-JD (universal-time->JD u)))

(defun hms->ra-rdn (hr min sec)
  "Return a double giving the number of radians right ascension corresponding to HR MIN, SEC time
hours, minutes, seconds.  Normally MIN is minute-of-hour and SEC is second-of-minute, but this isn't
essential; the returned value corresponds to the sum of the total times of HR, MIN and SEC."
  (* 2.0d0 pi (hms->day hr min sec)))

;;;; EARTH GEOMETRY

(defun xy-z->latdetic-altitude (xy z &key
                                (req 6378.135d0) ;Earth radius in km for WGS-72
                                (flattening (/ 298.26))
                                (lat-tol 1.0d-6) ;1 urad = about 7 m on the Earth's surface
                                (max-iter 20))
  "Return (values <geodetic latitude in radians> <altitude>) at the point XY, Z of the oblate body
with equatorial radius REQ and FLATTENING = 1 - (polar-radius/req).  Here XY is the radius in the xy
plane.  Solution is iterative.  If supplied, LAT-TOL is the tolerance on latitude in radians and
MAX-ITER is the maximum number of iterations allowed.  If these are not supplied, 'reasonable'
defaults will be used.  If convergence fails, an error is thrown."
  (let* ((esq (* flattening (- 2.0d0 flattening))) ;eccentricity squared
         (req-esq (* req esq)))
    (loop
       with latdetic = (atan (/ z xy)) ;initial estimate based on spherical shape
       for iter from 0
       for sin-lat = (sin latdetic)
       for c = (/ (sqrt (- 1.0d0 (* esq (* sin-lat sin-lat)))))
       for next-latdetic = (atan (/ (+ z (* (* req-esq c) sin-lat)) xy))
       do
       (cond ((< (abs (- next-latdetic latdetic)) lat-tol)
              (return (values next-latdetic (- (/ xy (cos next-latdetic)) (* req c)))))
             ((> iter max-iter)
              (error "convergence failed for xy = ~a, z = ~a" xy z))
             (t (setf latdetic next-latdetic))))))

(defun txyz->lla (time x y z &key
                  (f-time->phi #'phi-IAU-at-universal-time)
                  (req 6378.135d0) ;Earth radius in km for WGS-72
                  (flattening (/ 298.26))
                  (lat-tol 1.0d-6) ;1 urad = about 7 m on the Earth's surface
                  (max-iter 20))
  "F-TIME->PHI is a function that takes a TIME and returns the GMST or GHA in radians.  It
thus should be consistent with the scale used for TIME."
  (let ((xy (sqrt (+ (* x x) (* y y))))
        (gha (funcall f-time->phi time)))
    (multiple-value-bind (latdetic altitude)
        (xy-z->latdetic-altitude xy z
                                 :req req :flattening flattening :lat-tol lat-tol :max-iter max-iter)
      (let ((elong (rdn-mpi-ppi (- (atan y x) gha))))
        (values elong latdetic altitude)))))
