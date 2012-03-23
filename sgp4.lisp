;; -*- fill-column:102; comment-column:24; -*-

;; Copyright (c) Mayes Mullins, 2012
;; mmullins@mullinsenterprises.ca

(in-package :cl-sgp4)

;; (declaim (optimize (debug 3)))
(declaim (optimize (speed 3) (safety 0)))


;;;; UTILITIES

(defun make-v-el (n e I O w M)
  "Return a new orbital elements 6-vector containing the input parameters as values."
  (declare (type double-float n e I O w M))
  (let ((v-el (make-array 6 :element-type 'double-float)))
     (setf (aref v-el 0) n)
     (setf (aref v-el 1) e)
     (setf (aref v-el 2) I)
     (setf (aref v-el 3) O)
     (setf (aref v-el 4) w)
     (setf (aref v-el 5) M)
     v-el))

(defun make-v-el-fixing-negative-I (n e I O w M)
  "Return a new orbital elements 6-vector defining the same orbit as the discrete input values, but,
if the input set has a negative inclination, modify the element set to have a positive inclination."
  (if (>= I 0.0d0)
      (make-v-el n e I O w M)
      (make-v-el n e (- I) (+ O pi) (- w pi) M)))

(defun copy-v-el (v-el &rest i-v-pairs)
  "Return a copy of orbital elements 6-vector V-EL.  If more args are supplied, they should be pairs,
index and value to be substituted at that index in the copy.
Example: (copy-v-el #(0 1 2 3 4 5) 1 -1 4 -4) ==> #(0 -1 2 3 -4 5)"
  (let ((v-el-1 (make-array 6 :element-type 'double-float)))
    (loop for i from 0 below 6
       do (setf (aref v-el-1 i) (aref v-el i)))
    (loop for (i v) on i-v-pairs by #'cddr
       do (setf (aref v-el-1 i) v))
    v-el-1))

(declaim (inline make-v-el copy-v-el make-v-el-fixing-negative-I))


;;;; CONDITIONS

(define-condition altitude-lt-0-error (error)
  ((text :initarg :text :reader text)))

(define-condition Keplers-equation-convergence-error (error)
  ((text :initarg :text :reader text)))


;;;; SIDEREAL TIME

;; Two places below (in INTEGRATE-FROM-T50-SYNCHRONOUS-RESONANCE and
;; INTEGRATE-FROM-T50-HALF-DAY-RESONANCE) use PHI-SGP4-AT-T50 from :astrometry package.

;; Could replace these calls with calls to PHI-IAU-AT-T50 to use IAU sidereal time calculation,
;; perhaps for better consistency with some other calculation.  No significant difference in the
;; results though.


;;;; EARTH GRAVITY PROPERTIES

;; Each of these functions returns a 5-vector containing
;;    aE-km       ;equatorial radius of the Earth in km
;;    ke          ;in (Earth radii)^(3/2)/minute
;;    J2 J3 J4    ;un-normaliized zonal harmonic values
;; describing the Earth's gravity field.  The slightly different values represent different snapshots
;; of accuracy from different eras.  MAKE-V-GRAV-WG2-72 is the default used by Vallado2006 and by the
;; test cases checked here.

(defun make-v-grav-wgs-72-low-precision ()
  "Return Earth gravity parameters for low precision WGS 72 case."
  ;; STR3 values from Vallado2006, VI B and STR3
  (let ( ;; mu 398600.79964   ;not quite consistent with aE-km and ke
        (aE-km 6378.135d0)
        (J2 0.001082616d0)
        (J3 -0.00000253881d0)
        (J4 -0.00000165597d0)
        (ke 0.0743669161d0))
    (vector aE-km ke J2 J3 J4)))

(defun make-v-grav-wgs-72 ()
  "Return Earth gravity parameters for normal WGS 72 case."
  ;; Vallado2006, VI B
  (let* ((mu 398600.8d0)  ;mu - Earth gravitational parameter in km^3/s^2
         (aE-km 6378.135d0)
         (J2 0.001082616d0)
         (J3 -0.00000253881d0)
         (J4 -0.00000165597d0)
         (ke (/ 60.0d0 (sqrt (/ (* aE-km aE-km aE-km) mu)))))
    (vector aE-km ke J2 J3 J4)))

(defun make-v-grav-wgs-84 ()
  "Return Earth gravity parameters for WGS 84 case."
  ;; Vallado2006, VI B
  (let* ((mu 398600.5d0)  ;mu - Earth gravitational parameter in km^3/s^2
         (aE-km 6378.137d0)
         (J2 0.00108262998905d0)
         (J3 -0.00000253215306d0)
         (J4 -0.00000161098761d0)
         (ke (/ 60.0d0 (sqrt (/ (* aE-km aE-km aE-km) mu)))))
    (vector aE-km ke J2 J3 J4)))


;;;;; SGP4 INITIALIZATION

(defun n-Kozai->Brouwer (v-grav v-el0)
  "Return the Brouwer mean motion corresponding to the Kozai elements in input V-EL0.  Input V-GRAV is
the vector of gravity parameters in the form returned by one of the MAKE-V-GRAV-xxx functions.
Hoots2004 App B, section A."
  (let* ((ke (aref v-grav 1))
         (J2 (aref v-grav 2))
         (n0 (aref v-el0 0))
         (e0 (aref v-el0 1))
         (I0 (aref v-el0 2))
         (k2 (* 0.5d0 J2)) ;a_E in Earth radii is 1
         (cosI0 (cos I0))
         (cosI0to2 (* cosI0 cosI0))
         (e0to2 (* e0 e0))
         (delta-asq (/ (* 1.5d0 k2 (- (* 3.0d0 cosI0to2) 1.0d0))
                       (expt (- 1.0d0 e0to2) 1.5d0)))
         (a1 (expt (/ ke n0) (/ 2.0d0 3.0d0)))
         (delta1 (/ delta-asq (* a1 a1)))
         (delta1to2 (* delta1 delta1))
         (delta1to3 (* delta1 delta1to2))
         (a2 (* a1 (- 1.0d0 (+ (* (/ 3.0d0) delta1)
                               (+ delta1to2 (* (/ 134.0d0 81.0d0) delta1to3))))))
         (delta0 (/ delta-asq (* a2 a2))))
    (declare (type double-float
                   ke J2 n0 e0 I0 k2 cosI0 cosI0to2 e0to2 delta-asq a1 delta1 delta1to2 delta1to3 a2
                   delta0))
    ;; (format t "a1,a2,delta0,n0 ==> n0 = ~a ~a ~a ~a ~a~%" a1 a2 delta0 n0 (/ n0 (+ 1.0d0 delta0)))
    (/ n0 (+ 1.0d0 delta0))))

;; NOTE - all subsequent functions take Brouwer elements

(defun init-secular-atmospheric-drag (v-grav v-el0 Bstar)
  "Return initialization data for secular atmospheric drag calculation, (values q0 s v-Ci v-Di).  v-Ci
is [C1 C2 C3 C4 C5] and v-Di is [D2 D3 D4].  Input V-GRAV is the vector of gravity parameters in the
form returned by one of the MAKE-V-GRAV-xxx functions.  V-EL0 is the elements vector at epoch.  BSTAR
is the B-star parameter for the satellite.  See Hoots2004, App B, section A.1."
  (declare (type double-float Bstar))
  (let* ((aE-km (aref v-grav 0))
         (ke (aref v-grav 1))
         (J2 (aref v-grav 2))
         (J3 (aref v-grav 3))
         (n0 (aref v-el0 0))
         (e0 (aref v-el0 1))
         (I0 (aref v-el0 2))
         (w0 (aref v-el0 4))
         (k2 (* 0.5d0 J2)) ;{1 \over 2} J_2 a_E^2 but a_E in Earth radii is 1
         (A30 (- J3))      ;-J_3 a_E^3 but a_E in Earth radii is 1
         (sinI0 (sin I0))
         (a0 (expt (/ ke n0) (/ 2.0d0 3.0d0))) ;semi-major axis
         (rp (* a0 (- 1.0d0 e0)))              ;radius of perigee in Earth radii
         (perigee-km (* (- rp 1.0d0) aE-km))
         (q0 (+ 1.0d0 (/ 120.0d0 aE-km)))
         (s (cond ((<= perigee-km 98.0d0) (+ 1.0d0 (/ 20.0d0 aE-km)))
                  ((< perigee-km 156.0d0) (+ 1.0d0 (/ (- perigee-km 78.0d0) aE-km)))
                  (t (+ 1.0d0 (/ 78.0d0 aE-km)))))
         (theta (cos I0))
         (theta2 (* theta theta)) ;\theta^2
         (xi (/ (- a0 s)))
         (xi2 (* xi xi))                ;\xi^2
         (xi3 (* xi xi2))               ;\xi^3
         (xi4 (* xi2 xi2))              ;\xi^4
         (xi5 (* xi xi4))               ;\xi^5
         (beta0to2 (- 1.0d0 (* e0 e0))) ;\beta_0^2
         (eta (* a0 (* e0 xi))) ;\eta
         (eta2 (* eta eta))     ;\eta^2
         (eta3 (* eta eta2))    ;\eta^3
         (eta4 (* eta2 eta2))   ;\eta^4
         (q0ms (- q0 s))
         (q0ms2 (* q0ms q0ms))
         (q0ms4 (* q0ms2 q0ms2)) ;(q_0-s)^4
         (one-eta2 (abs (- 1.0d0 eta2))) ;without abs in onemetasqm72 calc, it becomes complex when perigee < 0;
                        ;also, have to use this for 1-\eta^2 in other places below to match Fortran program;
                        ;Fortran program does this as PSISQ  = DABS(1.0D0-ETASQ) in SGP4Init
         (Cfactor (* q0ms4 xi4 (expt one-eta2 -3.5d0)))
                        ;(q0-s)^4 \xi^4 (1-\eta^2)^{-{7 \over 2}} - part of several C formulas
                        ;(expt one-eta2 -3.5d0) = (1-\eta^2)^{-{7 \over 2}}  = Fortran program SGP4Init's COEF1
         (C2 (* Cfactor (* n0 (+ (* a0 (+ 1.0d0 (* 1.5d0 eta2) (* 4.0d0 e0 eta) (* e0 eta3)))
                                 (* 1.5d0 (/ (* k2 xi) one-eta2) (+ -0.5d0 (* 1.5d0 theta2))
                                    (+ 8.0d0 (+ (* 24.0d0 eta2) (* 3.0d0 eta4))))))))
         (C1 (* Bstar C2))
         (C1to2 (* C1 C1))       ;C_1^2
         (C1to3 (* C1 C1to2))    ;C_1^3
         (C1to4 (* C1to2 C1to2)) ;C_1^4
         (C3 (if (> e0 1.0d-4) ;Fortran program has this; also changes XMCOF which seems to be the
                        ;extra check I added to delta-M calc below
                 (/ (* q0ms4 xi5 A30 n0 sinI0) (* k2 e0))
                 0.0d0))
                        ;note -- aE = 1 in Earth radii
         (C4 (* 2.0d0 Cfactor n0 a0 beta0to2
                (- (+ (* 2.0d0 eta (+ 1.0d0 (* e0 eta))) (* 0.5d0 e0) (* 0.5d0 eta3))
                   (* (/ (* 2.0d0 k2 xi) (* a0 one-eta2))
                      (+ (* 3.0d0 (- 1.0d0 (* 3.0d0 theta2))
                            (- (+ 1.0d0 (* 1.5d0 eta2)) (* 2.0d0 e0 eta)
                               (* 0.5d0 e0 eta3)))
                         (* 0.75d0 (- 1.0d0 theta2)
                            (- (* 2.0d0 eta2) (* e0 eta) (* e0 eta3))
                            (cos (* 2.0d0 w0))))))))
         (C5 (* 2.0d0 Cfactor a0 beta0to2
                (+ 1.0d0 (* 2.75d0 eta (+ eta e0)) (* e0 eta3))))
         (D2 (* 4.0d0 a0 xi C1to2))
         (D3 (* (/ 4.0d0 3.0d0) a0 xi2 (+ (* 17.0d0 a0) s) C1to3))
         (D4 (* (/ 2.0d0 3.0d0) (* a0 a0) xi3 (+ (* 221.0d0 a0) (* 31.0d0 s)) C1to4)))
    (declare (type double-float
                   aE-km ke J2 J3 n0 e0 I0 w0 k2 A30 sinI0 a0 rp perigee-km q0 s theta theta2
                   xi xi2 xi3 xi4 xi5 beta0to2 eta eta2 eta3 eta4 q0ms q0ms2 q0ms4 one-eta2 Cfactor
                   C2 C1 C1to2 C1to3 C1to4 C3 C4 C5 D2 D3 D4))
    (values q0 s (vector C1 C2 C3 C4 C5) (vector D2 D3 D4))))

(defun init-Earth-zonal-harmonics (v-grav v-el0)
  "Return a 6-vector giving element derivatives due to secular Earth zonal harmonics.  Input V-GRAV is
the vector of gravity parameters in the form returned by one of the MAKE-V-GRAV-xxx functions.  V-EL0
is the elements vector at epoch.  See Hoots2004, App B, section A.2."
  (let* ((ke (aref v-grav 1))
         (J2 (aref v-grav 2))
         (J4 (aref v-grav 4))
         (n0 (aref v-el0 0))
         (e0 (aref v-el0 1))
         (I0 (aref v-el0 2))
         (k2 (* 0.5d0 J2))    ;{1 \over 2} J_2 a_E^2 but a_E in Earth radii is 1
         (k2to2 (* k2 k2))    ;k_2^2
         (k4 (* -0.375d0 J4)) ;-{3 \over 8} J_4 a_E^4 but a_E in Earth radii is 1
         (a0 (expt (/ ke n0) (/ 2.0d0 3.0d0))) ;semi-major axis
         (a0to2 (* a0 a0))                     ;a_0^2
         (theta (cos I0))
         (theta2 (* theta theta))           ;\theta^2
         (theta3 (* theta theta2))          ;\theta^3
         (theta4 (* theta theta3))          ;\theta^4
         (beta0 (sqrt (- 1.0d0 (* e0 e0)))) ;\beta_0
         (beta0to2 (* beta0 beta0))         ;\beta_0^2
         (beta0to3 (* beta0 beta0to2))      ;\beta_0^3
         ;; denominators of terms
         (x2a0to2beta0to3 (* 2.0d0 a0to2 beta0to3))
         (x2a0to2beta0to4 (* x2a0to2beta0to3 beta0))
         (x16a0to4beta0to7 (* x2a0to2beta0to3 x2a0to2beta0to3 4.0d0 beta0))
         (x16a0to4beta0to8 (* x16a0to4beta0to7 beta0))
         (x4a0to4beta0to8 (* 0.25d0 x16a0to4beta0to8))
         (xa0to2beta0to4 (* 0.5d0 x2a0to2beta0to4))
         (x2a0to4beta0to8 (* 0.5d0 x4a0to4beta0to8))
         (dMdt (* (+ (/ (* 3.0d0 k2 (- (* 3.0d0 theta2) 1.0d0))
                        x2a0to2beta0to3)
                     (/ (* 3.0d0 k2to2 (+ 13.0d0 (* -78.0d0 theta2)
                                          (* 137.0d0 theta4)))
                        x16a0to4beta0to7))
                  n0))
         (dwdt (* (+ (/ (* -3.0d0 k2 (- 1.0d0 (* 5.0d0 theta2))) x2a0to2beta0to4)
                     (/ (* 3.0d0 k2to2 (+ 7.0d0 (* -114.0d0 theta2) (* 395.0d0 theta4)))
                        x16a0to4beta0to8)
                     (/ (* 5.0d0 k4 (+ 3.0d0 (* -36.0d0 theta2) (* 49.0d0 theta4)))
                        x4a0to4beta0to8))
                  n0))
         (dOdt (* (+ (/ (* -3.0d0 k2 theta) xa0to2beta0to4)
                     (/ (* 3.0d0 k2to2 (- (* 4.0d0 theta) (* 19.0d0 theta3)))
                        x2a0to4beta0to8)
                     (/ (* 5.0d0 k4 theta (- 3.0d0 (* 7.0d0 theta2)))
                        x2a0to4beta0to8))
                  n0))
         ;; dndt = dedt = dIdt = 0
         )
    (declare (type double-float
                   ke J2 J4 n0 e0 I0 k2 k2to2 k4 a0 a0to2 theta theta2 theta3 theta4
                   beta0 beta0to2 beta0to3 x2a0to2beta0to3 x2a0to2beta0to4 x16a0to4beta0to7
                   x16a0to4beta0to8 x4a0to4beta0to8 xa0to2beta0to4 x2a0to4beta0to8 dMdt dwdt dOdt))
    (make-v-el 0.0d0 0.0d0 0.0d0 dOdt dwdt dMdt)))

;;;; INITIALIZATION FOR SECULAR AND LONG PERIOD COEFFICIENTS FOR LUNAR AND SOLAR GRAVITY

(defun v-deldt-3rd-body-perturbation (v-el0 nx sinIx cosIx Ox sinwx coswx Cx)
  "Return (values v-Xi v-Zi v-Zij v-deldt) where v-deldt is the 6-vector of element derviatives for
3rd-body (lunar, solar) perturbations.  Input V-EL0 is the satellite elements vector at epoch.  The
other inputs are parameters of the perturnbing body.  The following are its orbital parameters: NX is
its mean motion, SINIX and COSIX are the sin and cos of its inclination, OX is the RA of its ascending
node, SINWX and COSWX are the sin and cos of its argument of perigee.  CX is its perturbation
coefficient.  Note that we take sin/cos of Ix and wx instead of just Ix and wx because existing
program values are not quite consistent with 23.4441 deg, etc.  See Hoots2004, App A, section A."
  (declare (type double-float nx sinIx cosIx Ox sinwx coswx Cx))
  (let* ((n0 (aref v-el0 0))
         (e0 (aref v-el0 1))
         (I0 (aref v-el0 2))
         (O0 (aref v-el0 3))
         (w0 (aref v-el0 4))
         (cosI0 (cos I0))
         (sinI0 (sin I0))
         (cosw0 (cos w0))
         (sinw0 (sin w0))
         (O0mOx (- O0 Ox))
         (cosO0mOx (cos O0mOx))
         (sinO0mOx (sin O0mOx))

         (a1 (+ (* coswx cosO0mOx) (* sinwx cosIx sinO0mOx)))
         (a3 (+ (- (* sinwx cosO0mOx)) (* coswx cosIx sinO0mOx)))
         (a7 (+ (- (* coswx sinO0mOx)) (* sinwx cosIx cosO0mOx)))
         (a8 (* sinwx sinIx))
         (a9 (+ (* sinwx sinO0mOx) (* coswx cosIx cosO0mOx)))
         (a10 (* coswx sinIx))
         (a2 (+ (* cosI0 a7) (* sinI0 a8)))
         (a4 (+ (* cosI0 a9) (* sinI0 a10)))
         (a5 (+ (- (* sinI0 a7)) (* cosI0 a8)))
         (a6 (+ (- (* sinI0 a9)) (* cosI0 a10)))

         (X1 (+ (* a1 cosw0) (* a2 sinw0)))
         (X2 (+ (* a3 cosw0) (* a4 sinw0)))
         (X3 (+ (- (* a1 sinw0)) (* a2 cosw0)))
         (X4 (+ (- (* a3 sinw0)) (* a4 cosw0)))
         (X5 (* a5 sinw0))
         (X6 (* a6 sinw0))
         (X7 (* a5 cosw0))
         (X8 (* a6 cosw0))

         (e0sq (* e0 e0))
         (Z31 (- (* 12.0d0 X1 X1) (* 3.0d0 X3 X3)))
         (Z32 (- (* 24.0d0 X1 X2) (* 6.0d0 X3 X4)))
         (Z33 (- (* 12.0d0 X2 X2) (* 3.0d0 X4 X4)))
         (Z11 (+ (* -6.0d0 a1 a5) (* e0sq (- (* -24.0d0 X1 X7) (* 6.0d0 X3 X5)))))
         (Z12 (- (* -6.0d0 (+ (* a1 a6) (* a3 a5)))
                 (* e0sq (+ (* 24.0d0 (+ (* X2 X7) (* X1 X8)))
                            (* 6.0d0 (+ (* X3 X6) (* X4 X5)))))))
         (Z13 (+ (* -6.0d0 a3 a6) (* e0sq (- (* -24.0d0 X2 X8) (* 6.0d0 X4 X6)))))
         (Z21 (+ (* 6.0d0 a2 a5) (* e0sq (- (* 24.0d0 X1 X5) (* 6.0d0 X3 X7)))))
         (Z22 (+ (* 6.0d0 (+ (* a4 a5) (* a2 a6)))
                 (* e0sq (- (* 24.0d0 (+ (* X2 X5) (* X1 X6)))
                            (* 6.0d0 (+ (* X4 X7) (* X3 X8)))))))
         (Z23 (+ (* 6.0d0 a4 a6) (* e0sq (- (* 24.0d0 X2 X6) (* 6.0d0 X4 X8)))))

         (e0sqp1 (+ 1.0d0 e0sq))
         (Z1 (+  (* 6.0d0 (+ (* a1 a1) (* a2 a2))) (* e0sqp1 Z31)))
         (Z2 (+ (* 12.0d0 (+ (* a1 a3) (* a2 a4))) (* e0sqp1 Z32)))
         (Z3 (+  (* 6.0d0 (+ (* a3 a3) (* a4 a4))) (* e0sqp1 Z33)))
         (eta0 (sqrt (- 1.0d0 e0sq))) ;seems to be this, though I missed it in Hoots2004
         (I0lt3deg? (or (< I0 5.2359877d-2) (> I0 (- pi 5.2359877d-2))))
         (dndt 0.0d0)
         (dedt (* -15.0d0 Cx nx (/ (* e0 eta0) n0) (+ (* X1 X3) (* X2 X4))))
         (dIdt (* (/ (* (- Cx) nx) (* 2.0d0 (* n0 eta0))) (+ Z11 Z13)))
         (dOdt (if I0lt3deg?
                   0.0d0
                   (* (/ (* Cx nx) (* 2.0d0 (* n0 (* eta0 sinI0)))) (+ Z21 Z23))))
         ;; Fortran program SGP4UNIT.FOR line 692
         ;; BUT check might be supposed to be for I not I0 (?)
         (dwdt (if I0lt3deg?
                   (* (/ (* Cx (* nx eta0)) n0) (- (+ Z31 Z33) 6.0d0))
                   (- (* (/ (* Cx (* nx eta0)) n0) (- (+ Z31 Z33) 6.0d0))
                      (* dOdt cosI0))))
         (dMdt (* (/ (* (- Cx) nx) n0) (- (+ Z1 Z3) (+ 14.0d0 (* 6.0d0 e0sq)))))
         (v-deldt (make-v-el dndt dedt dIdt dOdt dwdt dMdt))
         (v-Xi (vector X1 X2 X3 X4 X5 X6 X7 X8))
         (v-Zi (vector Z1 Z2 Z3))
         (v-Zij (vector Z11 Z12 Z13 Z21 Z22 Z23 Z31 Z32 Z33)))
    (declare (type double-float
                   n0 e0 I0 O0 w0 cosI0 sinI0 cosw0 sinw0 O0mOx cosO0mOx sinO0mOx
                   a1 a3 a7 a8 a9 a10 a2 a4 a5 a6 X1 X2 X3 X4 X5 X6 X7 X8 e0sq
                   Z31 Z32 Z33 Z11 Z12 Z13 Z21 Z22 Z23 e0sqp1 Z1 Z2 Z3
                   eta0 dndt dedt dIdt dOdt dwdt dMdt))
    (values v-Xi v-Zi v-Zij v-deldt)))

(defun init-v-deldt-solar (v-el0)
  "Return (values v-Xi v-Zi v-Zij v-deldt) for solar perturbations.
Input T50 is the epoch time in days after 1950.0, V-EL0 is the elements 6-vector at this epoch.
See Hoots2004, App A, section A."
  (let ((n 1.19459d-5)            ;solar mean motion in rdn/min
        ;; I 0.4091767351668026 ;radians = 23.4441 deg
        ;; w 4.90822888509247 ;radians = 281.2208 deg
        ;; instead of these, the sin,cos values below are what are used
        ;; in the existing program
        (sinI 0.39785416d0)
        (cosI 0.91744867d0)
        (sinw -0.98088458d0)
        (cosw 0.1945905d0)
        (O 0.0d0)
        (C 2.9864797d-6)       ;C_s solar perurbation coefficient, radians/min
        )
    (declare (type double-float n sinI cosI sinw cosw O C))
    (v-deldt-3rd-body-perturbation v-el0 n sinI cosI O sinw cosw C)))

(defun init-v-deldt-lunar (t50 v-el0)
  "Return (values v-Xi v-Zi v-Zij v-deldt) for lunar perturbations.
Input T50 is the epoch time in days after 1950.0, V-EL0 is the elements 6-vector at this epoch.
See Hoots2004, App A, section A."
  ;; NOTE - description in Hoots2004, App A, section A implies lunar and solar perturbation values
  ;; calculated here are for epoch time; dscom uses time at epoch + tc; where dscom is called in
  ;; sgp4init; however, tc is set to 0.0; so seems this is independent of time depending only on epoch
  ;; time, and so can be called just once during initialization
  (declare (type double-float t50))
  (let* ((day      (+ t50 18261.5d0))
                        ;solar and lunar ephemeris uses Jan 0.5, 1900.0d0 = JD 2415020.0d0 epoch =
                        ;1899-12-31T12:00:00.000Z
         (Ome      (rdn-mpi-ppi (- 4.5236020d0 (* 9.2422029d-4 day))))
                        ;\Omega_m_epsilon
         (cosOme   (cos Ome)) ;cos(\Omega_m_\epsilon)
         (sinOme   (sin Ome)) ;sin(\Omega_m_\epsilon)
         (cosIm    (- 0.91375164d0 (* 0.03568096d0 cosOme)))
                        ;cos(I_m); here 0.913... = cos(\epsilon)
                        ;cos(I_m_\epsilon); 0.035... = sin(\epsilon)
                        ;sin(I_m_\epsilon)
         (sinIm    (sqrt (- 1.0d0 (* cosIm cosIm))))
                        ;sin(I_m) -- will sign be ok here ? should be given
                        ;range
         (sinOm    (/ (* 0.089683511d0 sinOme) sinIm))
                        ;sin(\Omega_m); here 0.0896... seems to be sin of
                        ;the moon's inclination relative to the ecliptic
                        ;I_m_\epsilon - Hoots2004 gives 5.145396374
                        ;deg --> sin of 0.089683448 - close enough?
         (cosOm    (sqrt (- 1.0d0 (* sinOm sinOm)))) ;cos(\Omega_m)
         (gamma    (+ 5.8351514d0 (* 0.0019443680d0 day)))
                        ;lunar longitude of perigee referred to the ecliptic
         (sinDelta (/ (* 0.39785416d0 sinOme) sinIm))
                        ;sin(\Delta); here 0.397... is sin(\epsilon)
         (cosDelta (+ (* cosOm cosOme) (* 0.91744867d0 sinOm sinOme)))
                        ;cos(\Delta); here 0.917... is cos(\epsilon)
         (Delta    (atan sinDelta cosDelta)) ;\Delta
         (n        1.583521770d-4)           ;lunar mean motion in rdn/min
         (O        (atan sinOm cosOm))       ;??? HOW TO KNOW SIGN OF cosOm ???
         (w        (- (+ gamma Delta) Ome))
                        ;\omega_m = \gamma - \Omega_m_\epsilon + \Delta =
                        ;G_0_m
         (sinw     (sin w))
         (cosw     (cos w))
         (C        4.7968065d-7)
                        ;C_m lunar perurbation coefficient, radians/min; this is program value, Hoots
                        ;has 4.796806521d-7
         )
    (declare (type double-float
                   day Ome cosOme sinOme cosIm sinIm sinOm cosOm gamma
                   sinDelta cosDelta Delta n O w sinw cosw C ))
    (v-deldt-3rd-body-perturbation v-el0 n sinIm cosIm O sinw cosw C)))

;;;; INITIALIZATION FOR RESONANCE EFFECTS OF EARTH GRAVITY

(defun init-v-deltai-synchronous-resonance (v-el0 ke)
  "Return a 3-vector containing [delta1 delta2 delta3], the delta values needed for synchronous
secular resonance calculations.  Input V-EL0 are the satellite elements at epoch, and KE is the k_e
gravity parameter of the Earth.  See Hoots2004, App A, section B."
  (declare (type double-float ke))
  (let* ((n0 (aref v-el0 0))
         (e0 (aref v-el0 1))
         (I0 (aref v-el0 2))
         (cosI0  (cos I0))
         (sinI0  (sin I0))
         (a0     (expt (/ ke n0) (/ 2.0d0 3.0d0)))
         (n0sq   (* n0 n0))
         (n2oa2  (/ n0sq (* a0 a0)))
         (n2oa3  (/ (* n0sq n0sq) (* ke ke))) ;n0^4/ke^2 = n0^2/a0^3
         (e0sq   (* e0 e0))
         (q22    1.7891679d-6)
         (q31    2.1460748d-6)
         (q33    2.2123015d-7)
         ;; Fijk are functions of epoch inclination
         (F220   (* 0.75d0 (* (+ 1.0d0 cosI0) (+ 1.0d0 cosI0))))
         (F311   (- (* 0.9375d0 (* sinI0 (* sinI0 (+ 1.0d0 (* 3.0d0 cosI0)))))
                    (* 0.75d0 (+ 1.0d0 cosI0))))
         (F330   (+ 1.0d0 cosI0))
         (F330   (* 1.875d0 (* F330 (* F330 F330))))
         ;; Gijk are functions of epoch eccentricity
         (G200   (+ 1.0d0 (* e0sq (+ -2.5d0 (* 0.8125d0 e0sq)))))
         (G310   (+ 1.0d0 (* 2.0d0 e0sq)))
         (G300   (+ 1.0d0 (* e0sq (+ -6.0d0 (* 6.60937d0 e0sq)))))
         ;; deltai
         (delta1 (* 3.0d0 n2oa3 (* F311 (* G310 q31))))
         (delta2 (* 6.0d0 n2oa2 (* F220 (* G200 q22))))
         (delta3 (* 9.0d0 n2oa3 (* F330 (* G300 q33)))))
    (declare (type double-float
                   n0 e0 I0 cosI0 sinI0 a0 n0sq n2oa2 n2oa3 e0sq q22 q31 q33
                   F220 F311 F330 F330 G200 G310 G300 delta1 delta2 delta3))
    (vector delta1 delta2 delta3)))

(defun init-v-Dlmpq-half-day-resonance (v-el0 ke)
  "Return a 10-vector containing the D_{lmpq} values needed for half-day secular resonance
calculations.  These are in the order
  [D2201 D2211 D3210 D3222 D4410 D4422 D5220 D5232 D5421 D5433]
Input V-EL0 are the satellite elements at epoch, and KE is the k_e gravity parameter of the Earth.
See Hoots2004, App A, section B."
  (declare (type double-float ke))
  (let* ((n0 (aref v-el0 0))
         (e0 (aref v-el0 1))
         (I0 (aref v-el0 2))
         ;; Fijk are functions of epoch inclination
         (cosI0   (cos I0))
         (sinI0   (sin I0))
         (cosI0sq (* cosI0 cosI0))
         (sinI0sq (- 1.0d0 cosI0sq))
         (F220    (* 0.75d0 (+ 1.0d0 (* 2.0d0 cosI0) cosI0sq)))
         (F221    (* 1.5d0 sinI0sq))
         (F321    (* 1.875d0 (* sinI0 (- 1.0d0 (* 2.0d0 cosI0) (* 3.0d0 cosI0sq)))))
         (F322    (* -1.875d0 (* sinI0 (- (+ 1.0d0 (* 2.0d0 cosI0)) (* 3.0d0 cosI0sq)))))
         (F441    (* 35.0d0 (* sinI0sq F220)))
         (F442    (* 39.375d0 (* sinI0sq sinI0sq)))
         (F522    (* 9.84375d0
                     (* sinI0
                        (+ (* sinI0sq (- 1.0d0 (* 2.0d0 cosI0) (* 5.0d0 cosI0sq)))
                           (*
                            0.33333333d0 ;; (/ 3.0d0)
                            (+ -2.0d0 (* 4.0d0 cosI0) (* 6.0d0 cosI0sq)))))))
         (F523    (* sinI0
                     (+ (* 4.92187512d0 (* sinI0sq (+ (- -2.0d0 (* 4.0d0 cosI0)) (* 10.0d0 cosI0sq))))
                        (* 6.56250012d0 (- (+ 1.0d0 (* 2.0d0 cosI0)) (* 3.0d0 cosI0sq))))))
         (F542    (* 29.53125d0
                     (* sinI0 (+ (- 2.0d0 (* 8.0d0 cosI0))
                                 (* cosI0sq (+ -12.0d0 (* 8.0d0 cosI0) (* 10.0d0 cosI0sq)))))))
         (F543    (* 29.53125d0
                     (* sinI0 (+ (- -2.0d0 (* 8.0d0 cosI0))
                                 (* cosI0sq (- (+ 12.0d0 (* 8.0d0 cosI0)) (* 10.0d0 cosI0sq)))))))
         ;; Gijk are functions of epoch eccentricity
         (e0sq    (* e0 e0))
         (e0cu    (* e0 e0sq))
         (G201    (- -0.306d0 (* (- e0 0.64d0) 0.44d0)))
         ;; would be more efficient to gang up these tests instead of having one test for each Gxyz
         (G211    (if (<= e0 0.65d0)
                      (+ (- 3.616d0 (* 13.247d0 e0)) (* 16.29d0 e0sq))
                      (+ (- (+ -72.099d0 (* 331.819d0 e0)) (* 508.738d0 e0sq))
                         (* 266.724d0 e0cu))))
         (G310    (if (<= e0 0.65d0)
                      (+ (- (+ -19.302d0 (* 117.39d0 e0)) (* 228.419d0 e0sq))
                         (* 156.591d0 e0cu))
                      (+ (- (+ -346.844d0 (* 1582.851d0 e0)) (* 2415.925d0 e0sq))
                         (* 1246.113d0 e0cu))))
         (G322    (if (<= e0 0.65d0)
                      (+ (- (+ -18.9068d0 (* 109.7927d0 e0)) (* 214.6334d0 e0sq))
                         (* 146.5816d0 e0cu))
                      (+ (- (+ -342.585d0 (* 1554.908d0 e0)) (* 2366.899d0 e0sq))
                         (* 1215.972d0 e0cu))))
         (G410    (if (<= e0 0.65d0)
                      (+ (- (+ -41.122d0 (* 242.694d0 e0)) (* 471.094d0 e0sq))
                         (* 313.953d0 e0cu))
                      (+ (- (+ -1052.797d0 (* 4758.686d0 e0)) (* 7193.992d0 e0sq))
                         (* 3651.957d0 e0cu))))
         (G422    (if (<= e0 0.65d0)
                      (+ (- (+ -146.407d0 (* 841.88d0 e0)) (* 1629.014d0 e0sq))
                         (* 1083.435d0 e0cu))
                      (+ (- (+ -3581.69d0 (* 16178.11d0 e0)) (* 24462.77d0 e0sq))
                         (* 12422.52d0 e0cu))))
         (G520    (if (<= e0 0.65d0)
                      (+ (- (+ -532.114d0 (* 3017.977d0 e0)) (* 5740.032d0 e0sq))
                         (* 3708.276d0 e0cu))
                      (if (> e0 0.715d0) ;Hoots has >=, existing program has >
                          (+ (- (+ -5149.66d0 (* 29936.92d0 e0)) (* 54087.36d0 e0sq))
                             (* 31324.56d0 e0cu))
                          (+ (- 1464.74d0 (* 4664.75d0 e0)) (* 3763.64d0 e0sq)))))
         (G533    (if (< e0 0.7d0)
                      (+ (- (+ -919.2277d0 (* 4988.61d0 e0)) (* 9064.77d0 e0sq))
                         (* 5542.21d0 e0cu))
                      (+ (- (+ -37995.78d0 (* 161616.52d0 e0)) (* 229838.2d0 e0sq))
                         (* 109377.94d0 e0cu))))
         (G521    (if (< e0 0.7d0)
                      (+ (- (+ -822.71072d0 (* 4568.6173d0 e0)) (* 8491.4146d0 e0sq))
                         (* 5337.524d0 e0cu))
                      (+ (- (+ -51752.104d0 (* 218913.95d0 e0)) (* 309468.16d0 e0sq))
                         (* 146349.42d0 e0cu))))
         (G532    (if (< e0 0.7d0)
                      (+ (- (+ -853.666d0 (* 4690.25d0 e0)) (* 8624.77d0 e0sq))
                         (* 5341.4d0 e0cu))
                      (+ (- (+ -40023.88d0 (* 170470.89d0 e0)) (* 242699.48d0 e0sq))
                         (* 115605.82d0 e0cu))))
         ;; Dlmpq
         (root22  1.7891679d-6)        ;\sqrt( C_{22}^2 + S_{22}^2 ), etc
         (root44  7.3636953d-9)
         (root54  2.1765803d-9)
         (root32  3.7393792d-7)
         (root52  1.1428639d-7)
         (a-recip (expt (/ n0 ke) (/ 2.0d0 3.0d0)))
         (n0sq    (* n0 n0))
         (coeff1  (* 3.0d0 (* n0sq (* a-recip a-recip))))
         (coeff   (* coeff1 root22))
         (D2201   (* coeff (* F220 G201)))
         (D2211   (* coeff (* F221 G211)))
         (coeff1  (* coeff1 a-recip))
         (coeff   (* coeff1 root32))
         (D3210   (* coeff (* F321 G310)))
         (D3222   (* coeff (* F322 G322)))
         (coeff1  (* coeff1 a-recip))
         (coeff   (* 2.0d0 (* coeff1 root44)))
         (D4410   (* coeff (* F441 G410)))
         (D4422   (* coeff (* F442 G422)))
         (coeff1  (* coeff1 a-recip))
         (coeff   (* coeff1 root52))
         (D5220   (* coeff (* F522 G520)))
         (D5232   (* coeff (* F523 G532)))
         (coeff   (* 2.0d0 (* coeff1 root54)))
         (D5421   (* coeff (* F542 G521)))
         (D5433   (* coeff (* F543 G533))))
    (declare (type double-float
                   n0 e0 I0 cosI0 sinI0 cosI0sq sinI0sq
                   F220 F221 F321 F322 F441 F442 F522 F523 F542 F543 e0sq e0cu
                   G201 G211 G310 G322 G410 G422 G520 G533 G521 G532
                   root22 root44 root54 root32 root52 a-recip n0sq coeff1 coeff
                   D2201 D2211 D3210 D3222 D4410 D4422 D5220 D5232 D5421 D5433 ))
    (vector D2201 D2211 D3210 D3222 D4410 D4422 D5220 D5232 D5421 D5433)))


;;;; USER LEVEL INITIALIZATION DEFSTRUCT AND FUNCTION

(defstruct sgp4data

  ;; initialization data derived from epoch time and tledata; user can treat
  ;; this as a black box, produced by function init-sgp4 and used by
  ;; function propagate-sgp4

  v-grav                        ;[aE-km ke J2 J3 J4]
  t50                           ;epoch in days after 1950.0
  v-el0                         ;[n e I O w M] elements at epoch
                                ;Brouwer ("un-Kozai'd") values
  Bstar                         ;drag parameter B-star
  q0                            ;drag parameter q_0
  s                             ;drag parameter s
  v-Ci                          ;[C1 C2 C3 C4 C5]   ???? C2 not used after making C1
  v-Di                          ;[D2 D3 D4]
  v-deldt-Earth-zonal-harmonics ;element derivatives for Earth harmonics
  deep-space?                   ;if true, deep-space calculations are done
  resonance                     ;:synchronous or :half-day or nil
  drop-terms?                   ;if true, some higher order terms in drag
                                ;calculations are dropped
                                ;; fields v-Xi-lunar ... v-deldt-solar are set if deep-space? is true
                                ;; are nil otherwise
  v-Xi-lunar                    ;lunar [X1 X2 X3 X4 X5 X6 X7 X8]
  v-Zi-lunar                    ;lunar [Z1 Z2 Z3]
  v-Zij-lunar                   ;lunar [Z11 Z12 Z13 Z21 Z22 Z23 Z31 Z32 Z33]
  v-deldt-lunar                 ;element derivatives for lunar perturbations
  v-Xi-solar                    ;solar [X1 X2 X3 X4 X5 X6 X7 X8]
  v-Zi-solar                    ;solar [Z1 Z2 Z3]
  v-Zij-solar                   ;solar [Z11 Z12 Z13 Z21 Z22 Z23 Z31 Z32 Z33]
  v-deldt-solar                 ;element derivatives for solar perturbations
  v-deltai                      ;[delta1 delta2 delta3]
                                ;set if resonance == :synchronous, nil otherwise
  v-Dlmpq                       ;[D2201 D2211 D3210 D3222 D4410 D4422
                                ; D5220 D5232 D5421 D5433]
                                ;set if resonance == :half-day, nil otherwise
  )

(defun init-sgp4 (v-grav t50 v-el0 Bstar)
  "Do initialization calcs and return an SGP4DATA initialization record for an SGP4 satellite.  Input
V-GRAV is the vector of gravity parameters in the form returned by one of the MAKE-V-GRAV-xxx
functions.  T50 is the epoch time in days after 1950.0.  V-EL0 is the elements vector at epoch.
BSTAR is the B-star parameter for the satellite."
  (declare (type double-float t50 Bstar))
  (let* ((aE-km (aref v-grav 0))
         (ke (aref v-grav 1))
         (e0 (aref v-el0 1))
         (n0B (n-Kozai->Brouwer v-grav v-el0))
         (v-el0B (copy-v-el v-el0 0 n0B)) ;"Brouwer"/"un-Kozai'd" elements at epoch
         (v-deldt-Earth-zonal-harmonics (init-Earth-zonal-harmonics v-grav v-el0B))
         ;; categorize satellite
         (period-min (/ (* 2.0d0 pi) n0B)) ;n0B is in radians/min
         (deep-space? (>= period-min 225.0d0))
         (resonance (cond ((and (>= period-min 1200.0d0) (<= period-min 1800.0d0))
                        ;dsinit uses ((nm < 0.0052359877) && (nm >
                        ;0.0034906585)) where nm = n0B here; Hoots2004
                        ;calls it a "closed interval", so I've used <=/>=
                           :synchronous)
                          ((and (>= period-min 680.0d0) (<= period-min 760.0d0)
                                (>= e0 0.5d0))
                        ;dsinit uses ((nm >= 8.26d-3) && (nm <= 9.24d-3) &&
                        ;(em >= 0.5)); eccentricity check is not in Hoots
                           :half-day)
                          (t nil)))
         (a0 (expt (/ ke n0B) (/ 2.0d0 3.0d0))) ;semi-major axis
         (rp (* a0 (- 1.0d0 e0)))            ;radius of perigee in Earth radii
         (perigee-km (* (- rp 1.0d0) aE-km))
         (drop-terms? (or (< perigee-km 220.0d0) deep-space?))
         ;; make values used only if there is a resonance
         (v-deltai (when (eq resonance :synchronous)
                     (init-v-deltai-synchronous-resonance v-el0B ke)))
         (v-Dlmpq (when (eq resonance :half-day)
                    (init-v-Dlmpq-half-day-resonance v-el0B ke))))
    (declare (type double-float aE-km ke e0 n0B period-min a0 rp perigee-km))
    (multiple-value-bind (q0 s v-Ci v-Di) (init-secular-atmospheric-drag v-grav v-el0B Bstar)
      ;; make values used only if deep-space?
      (multiple-value-bind (v-Xi-solar v-Zi-solar v-Zij-solar v-deldt-solar)
          (if deep-space?
              (init-v-deldt-solar v-el0B)
              (values nil nil nil nil))
        (multiple-value-bind (v-Xi-lunar v-Zi-lunar v-Zij-lunar v-deldt-lunar)
            (if deep-space?
                (init-v-deldt-lunar t50 v-el0B)
                (values nil nil nil nil))
          (make-sgp4data
           :v-grav v-grav :t50 t50 :v-el0 v-el0B :Bstar Bstar :q0 q0 :s s :v-Ci v-Ci :v-Di v-Di
           :v-deldt-Earth-zonal-harmonics v-deldt-Earth-zonal-harmonics
           :deep-space? deep-space? :resonance resonance :drop-terms? drop-terms?
           :v-Xi-lunar v-Xi-lunar :v-Zi-lunar v-Zi-lunar :v-Zij-lunar v-Zij-lunar :v-deldt-lunar v-deldt-lunar
           :v-Xi-solar v-Xi-solar :v-Zi-solar v-Zi-solar :v-Zij-solar v-Zij-solar :v-deldt-solar v-deldt-solar
           :v-deltai v-deltai :v-Dlmpq v-Dlmpq))))))


;;;;; SGP4 PROPAGATION

;;;; SECULAR UPDATE FOR EARTH ZONAL GRAVITY AND PARTIAL ATMOSPHERIC DRAG

(defun update-v-el-for-Earth-zonal-gravity-and-partial-drag
    (v-el0 v-grav Bstar q0 s v-Ci v-deldt-Earth-zonal-harmonics drop-terms? time)
  "Return an updated element set 6-vector accounting for Earth zonal gravity and partial drag effects.
Inputs are components of an SGP4DATA record plus TIME which is the time in minutes after epoch.  See
Hoots2004, App B, section B.1."
  (declare (type double-float Bstar q0 s time))
  (let* (;; (aE-km (aref v-grav 0))
         (ke (aref v-grav 1))
         (J2 (aref v-grav 2))
         (n0 (aref v-el0 0))
         (e0 (aref v-el0 1))
         (I0 (aref v-el0 2))
         (O0 (aref v-el0 3))
         (w0 (aref v-el0 4))
         (M0 (aref v-el0 5))
         (C1 (aref v-Ci 0))
         (C3 (aref v-Ci 2))
         (dO0dt (aref v-deldt-Earth-zonal-harmonics 3))
         (dw0dt (aref v-deldt-Earth-zonal-harmonics 4))
         (dM0dt (aref v-deldt-Earth-zonal-harmonics 5))
         ;; could generate the following constants once and save in sgp4data; they're all in
         ;; init-secular-atmospheric-drag; the ones needed here are k2, not a0, the ones with TeX
         ;; formula comments; if I'm going to not save them, but re-calculate them here, why bother
         ;; with q0 and s then ?  they're just as easy to calculate;
         (e0m (max e0 1.0d-6)) ;fix for near zero eccentricity
         (aE 1.0d0)            ;radius of the Earth in Earth radii units
         (k2 (* 0.5d0 J2))     ;{1 \over 2} J_2 a_E^2 but a_E in Earth radii is 1
         (a0 (expt (/ ke n0) (/ 2.0d0 3.0d0))) ;semi-major axis
         (a0to2 (* a0 a0))                ;a_0^2
         (theta (cos I0))                 ;\theta
         (beta0to2 (- 1.0d0 (* e0m e0m)))     ;\beta_0^2
         (q0ms (- q0 s))
         (q0ms2 (* q0ms q0ms))
         (q0ms4 (* q0ms2 q0ms2))       ;(q_0-s)^4
         (xi (/ (- a0 s)))
         (xi2 (* xi xi))
         (xi4 (* xi2 xi2))             ;\xi^4
         (eta (* a0 (* e0m xi)))        ;\eta
         (Mdf (+ M0 (* n0 time) (* dM0dt time)))
         (wdf (+ w0 (* dw0dt time)))
         (Odf (+ O0 (* dO0dt time)))
         (delta-w (if drop-terms?
                      0.0d0
                      (* Bstar C3 (cos w0) time)))
         (delta-M (if (or (< e0m 1.0d-4) drop-terms?) ;e0m check is in Fortran program, not Hoots
                      0.0d0
                      (* (/ -2.0d0 3.0d0) q0ms4 Bstar xi4
                         (/ aE (* e0m eta))
                         (- (expt (+ 1.0d0 (* eta (cos Mdf))) 3)
                            (expt (+ 1.0d0 (* eta (cos M0))) 3)))))
         (M (+ Mdf delta-w delta-M))
         (w (- wdf (+ delta-w delta-M)))
         (O (- Odf (* 10.5d0 (/ (* n0 k2 theta)
                                  (* a0to2 beta0to2)) C1 (* time time)))))
    (declare (type double-float
                   ke J2 n0 e0 I0 O0 w0 M0 C1 C3 dO0dt dw0dt dM0dt e0m aE k2 a0 a0to2
                   theta beta0to2 q0ms q0ms2 q0ms4 xi xi2 xi4 eta Mdf wdf Odf delta-w delta-M M w O ))
    (make-v-el n0 e0m I0 O w M)))


;;;; SECULAR UPDATES FOR EFFECTS OF LUNAR AND SOLAR GRAVITY

(defun update-v-el-for-secular-lunar-solar-gravity (v-el0 v-deldt-lunar v-deldt-solar time)
  "Return a new element set updating V-EL0 for lunar and solar secular rates in V-DELDT-LUNAR and
V-DELDT-SOLAR and time TIME minutes after epoch.  See Hoots2004, App B, section B.2."
  (declare (type double-float time))
  (let ((v-el (make-array 6 :element-type 'double-float)))
    (loop for i below 6
       do (setf (aref v-el i) (+ (aref v-el0 i) (* time (+ (aref v-deldt-lunar i)
                                                           (aref v-deldt-solar i))))))
    v-el))


;;;; SECULAR UPDATES FOR RESONANCE EFFECTS OF EARTH GRAVITY

(defun derivs-synchronous-resonance
    (v-deltai     ;[\delta_1 \delta_2 \delta_3]
     lambda0dot   ;\lambda_0\dot
     n            ;mean motion
     lambda       ;\lambda = "resonance variable" in Euler-Maclaurin integration
     )
  "Return (value dndt dldt d2ndt2 d2ldt2 giving the first and second derivatives of N (mean motion)
and LAMBDA for secular resonance effects of Earth gravity on near synchronous satellites.  See
Hoots2004 2004, App A, section D."
  (declare (type double-float lambda0dot n lambda))
  (let* ((lambda31 0.13130908d0) ;\lambda_{31} App A B p.4
         (lambda22 2.88431980d0) ;\lambda_{22}
         (lambda33 0.37448087d0) ;\lambda_{33}
         (delta1 (aref v-deltai 0))
         (delta2 (aref v-deltai 1))
         (delta3 (aref v-deltai 2))
         (dndt     (+ (* delta1 (sin (- lambda lambda31)))
                      (* delta2 (sin (* 2.0d0 (- lambda lambda22))))
                      (* delta3 (sin (* 3.0d0 (- lambda lambda33))))))
         (dldt     (+ n lambda0dot))
         (d2ndt2   (* dldt (+
                              (* delta1 (cos (- lambda lambda31)))
                              (* 2.0d0 (* delta2
                                          (cos (* 2.0d0 (- lambda lambda22)))))
                              (* 3.0d0 (* delta3
                                          (cos (* 3.0d0 (- lambda lambda33))))))))
         (d2ldt2 dndt))
    (declare (type double-float
                   lambda31 lambda22 lambda33 delta1 delta2 delta3 dndt dldt d2ndt2 d2ldt2))
    (values dndt dldt d2ndt2 d2ldt2)))

(defun derivs-half-day-resonance
    (v-Dlmpq        ;[D2201 D2211 D3210 D3222 D4410 D4422 D5220 D5232 D5421 D5433]
                        ;these are D_{lmpq}, i.e., D2201 = D_{2201}, etc
     lambda0dot                           ;\lambda_0\dot
     omega          ;\omega = argument of perigee at E-M step i
                        ;"secularly updated argument of perigee
                        ;\omega_0 + \omega_0\dot \Delta time"
     n                                    ;mean motion
     lambda         ;\lambda = "resonance variable" in Euler-Maclaurin integration
     )
  "Return (values dndt dldt d2ndt2 d2ldt2) giving the first and second derivatives of N (mean motion)
and LAMBDA for secular resonance effects of Earth gravity on near half-day satellites.  See Hoots2004
2004, App A, section D."
  (declare (type double-float lambda0dot omega n lambda))
  (let* ((g22     5.7686396d0)             ;G_{22} App A D p.6
         (g32     0.95240898d0)    ;G_{32}
         (g44     1.8014998d0)     ;G_{44}
         (g52     1.0508330d0)     ;G_{52}
         (g54     4.4108898d0)     ;G_{54}
         (D2201 (aref v-Dlmpq 0))
         (D2211 (aref v-Dlmpq 1))
         (D3210 (aref v-Dlmpq 2))
         (D3222 (aref v-Dlmpq 3))
         (D4410 (aref v-Dlmpq 4))
         (D4422 (aref v-Dlmpq 5))
         (D5220 (aref v-Dlmpq 6))
         (D5232 (aref v-Dlmpq 7))
         (D5421 (aref v-Dlmpq 8))
         (D5433 (aref v-Dlmpq 9))
         (omega2  (+ omega omega)) ;2\omega
                        ;(not \omega^2 as in most other uses of '2' suffix)
         (lambda2 (+ lambda lambda)) ;2\lambda (ditto)
         (omega+lambda (+ omega lambda))
         (omega2+lambda (+ omega2 lambda))
         (omega+lambda2 (+ omega lambda2))
         (omega2+lambda2 (+ omega2 lambda2))
         (dndt    (+ (* D2201 (sin (- omega2+lambda g22)))
                     (* D2211 (sin (- lambda g22)))
                     (* D3210 (sin (- omega+lambda g32)))
                     (* D3222 (sin (- lambda omega g32)))
                     (* D4410 (sin (- omega2+lambda2 g44)))
                     (* D4422 (sin (- lambda2 g44)))
                     (* D5220 (sin (- omega+lambda g52)))
                     (* D5232 (sin (- lambda omega g52)))
                     (* D5421 (sin (- omega+lambda2 g54)))
                     (* D5433 (sin (- lambda2 omega g54)))))
         (dldt    (+ n lambda0dot))
         (d2ndt2  (* dldt (+ (* D2201 (cos (- omega2+lambda g22)))
                             (* D2211 (cos (- lambda g22)))
                             (* D3210 (cos (- omega+lambda g32)))
                             (* D3222 (cos (- (+ (- omega) lambda) g32)))
                             (* D5220 (cos (- omega+lambda g52)))
                             (* D5232 (cos (- (+ (- omega) lambda) g52)))
                             (* 2.0d0 (+ (* D4410 (cos (- omega2+lambda2 g44)))
                                         (* D4422 (cos (- lambda2 g44)))
                                         (* D5421 (cos (- omega+lambda2 g54)))
                                         (* D5433 (cos (- (+ (- omega) lambda2) g54))))))))
         (d2ldt2 dndt))
    (declare (type double-float
                   g22 g32 g44 g52 g54
                   D2201 D2211 D3210 D3222 D4410 D4422 D5220 D5232 D5421 D5433
                   omega2 lambda2 omega+lambda omega2+lambda omega+lambda2 omega2+lambda2
                   dndt dldt d2ndt2 d2ldt2))
    (values dndt dldt d2ndt2 d2ldt2)))

(defun integrate-synchronous-resonance
    (v-deltai           ;[\delta_1 \delta_2 \delta_3]
     lambda0dot         ;\lambda_0\dot
     time-prev n-prev l-prev
     dndt-prev dldt-prev d2ndt2-prev d2ldt2-prev
     time-final               ;final time for integration
     step               ;magnitude of step size to use in integration
     )
  "Do Euler-Maclaurin integration for Earth gravity synchronous resonance effects.
Return 2 values:
 - n at time TIME minutes after epoch
 - l (lambda) at time TIME minutes after epoch"
  ;; note -- not trying to make the most accurate integration; rather trying
  ;; to match what's done in current SGP4 model code; that uses a fixed step
  ;; length and 2nd order extrapolation from the previous step to the time
  ;; required
  (declare (type double-float lambda0dot
                 time-prev n-prev l-prev
                 dndt-prev dldt-prev d2ndt2-prev d2ldt2-prev
                 time-final step))
  (let* ((step (if (>= time-final time-prev) (abs step) (- (abs step)))))
    (declare (type double-float step))
    (loop
       with time of-type double-float = time-prev
       with n of-type double-float = n-prev
       with l of-type double-float = l-prev
       with dndt of-type double-float = dndt-prev
       with dldt of-type double-float = dldt-prev
       with d2ndt2 of-type double-float = d2ndt2-prev
       with d2ldt2 of-type double-float = d2ldt2-prev
       ;; note - these are altered each pass through the loop
       do (cond ((or (and (> step 0.0d0) (< time-final (+ time step)))
                     (and (< step 0.0d0) (> time-final (+ time step))))
                 (let ((dt (- time-final time)))
                   (return
                     (values
                      (+ n (* dt (+ dndt (* 0.5d0 dt d2ndt2))))
                      (+ l (* dt (+ dldt (* 0.5d0 dt d2ldt2))))))))
                (t      ;do a step
                 (incf time step)
                 (incf n (* step (+ dndt (* 0.5d0 step d2ndt2))))
                 (incf l (* step (+ dldt (* 0.5d0 step d2ldt2))))
                 (multiple-value-setq (dndt dldt d2ndt2 d2ldt2)
                   (derivs-synchronous-resonance v-deltai lambda0dot n l)))))))

(defun integrate-half-day-resonance
    (v-Dlmpq      ;[D2201 D2211 D3210 D3222 D4410 D4422 D5220 D5232 D5421 D5433]
     w0           ;\omega_0 = \omega at epoch, time_0
     dwdt0        ;d \omega / d time at epoch, time_0
     lambda0dot   ;\lambda_0\dot
     time-prev n-prev l-prev
     dndt-prev dldt-prev d2ndt2-prev d2ldt2-prev
     time-final   ;final time for integration (actually time-t50 = \Delta time)
     step         ;magnitude of step size to use in Euler-Maclaurin integration
     )
  "DoEuler-Maclaurin integration for Earth gravity half-day resonance effects.
Return 2 values:
 - n at time TIME minutes after epoch
 - l (lambda) at time TIME minutes after epoch "
  ;; note -- not trying to make the most accurate integration; rather trying
  ;; to match what's done in current SGP4 model code; that uses a fixed step
  ;; length and 2nd order extrapolation from the previous step to the time
  ;; required
  (declare (type double-float w0 dwdt0 lambda0dot
                 time-prev n-prev l-prev
                 dndt-prev dldt-prev d2ndt2-prev d2ldt2-prev
                 time-final step))
  (let* ((step (if (>= time-final time-prev) (abs step) (- (abs step)))))
    (declare (type double-float step))
    (loop
       with time of-type double-float = time-prev
       with n of-type double-float = n-prev
       with l of-type double-float = l-prev
       with dndt of-type double-float = dndt-prev
       with dldt of-type double-float = dldt-prev
       with d2ndt2 of-type double-float = d2ndt2-prev
       with d2ldt2 of-type double-float = d2ldt2-prev
       ;; note - these are altered each pass through the loop
       do (cond ((or (and (> step 0.0d0) (< time-final (+ time step)))
                     (and (< step 0.0d0) (> time-final (+ time step))))
                 (let ((dt (- time-final time))) ;remaining time after prev step
                   (return
                     (values (+ n (* dt (+ dndt (* 0.5d0 dt d2ndt2))))
                             (+ l (* dt (+ dldt (* 0.5d0 dt d2ldt2))))))))
                (t      ;do a step
                 (incf time step)
                 (incf n (* step (+ dndt (* 0.5d0 step d2ndt2))))
                 (incf l (* step (+ dldt (* 0.5d0 step d2ldt2))))
                 (multiple-value-setq (dndt dldt d2ndt2 d2ldt2)
                   (derivs-half-day-resonance v-Dlmpq lambda0dot
                                              (+ w0 (* dwdt0 time)) ;\omega_i = "secularly updated
                                              n l)))))))

(defun integrate-from-t50-synchronous-resonance
    (t50 v-el0
     v-deldt-Earth-zonal-harmonics
     v-deldt-lunar v-deldt-solar v-deltai
     time v-el)
  "Do Euler-Maclaurin integration for Earth gravity synchronous resonance effects starting from epoch
time and going to time TIME minutes after epoch.  V-EL are the elements at time TIME after accounting
for Earth zonal gravity, partial drag and lunar/solar gravity secular effects.  Return (values n M)."
  (declare (type double-float t50 time))
  (let* ((n0 (aref v-el0 0))
         (O0 (aref v-el0 3))
         (w0 (aref v-el0 4))
         (M0 (aref v-el0 5))
         (O (aref v-el 3))
         (w (aref v-el 4))
         (v-deldt-EZ v-deldt-Earth-zonal-harmonics)
         (dOdt-EZ (aref v-deldt-EZ 3))
         (dwdt-EZ (aref v-deldt-EZ 4))
         (dMdt-EZ (aref v-deldt-EZ 5))
         (dOdt-LS (+ (aref v-deldt-lunar 3) (aref v-deldt-solar 3)))
         (dwdt-LS (+ (aref v-deldt-lunar 4) (aref v-deldt-solar 4)))
         (dMdt-LS (+ (aref v-deldt-lunar 5) (aref v-deldt-solar 5)))
         (theta0 (phi-sgp4-at-t50 t50))
         (dthetadt 4.37526908801129966d-3) ;gmst/sidereal rotation rate in radians/minute
         (thetat (+ theta0 (* dthetadt time))) ;should be close enough
         (lambda0 (rdn-mpi-ppi (+ M0 (+ w0 (- O0 theta0)))))
         (lambda0dot (- (+ dMdt-EZ dMdt-LS (+ dOdt-EZ dOdt-LS)
                           (+ dwdt-EZ dwdt-LS)) dthetadt))
         (step 720.0d0))
    (declare (type double-float
                   n0 O0 w0 M0 O w dOdt-EZ dwdt-EZ dMdt-EZ dOdt-LS dwdt-LS dMdt-LS
                   theta0 dthetadt thetat lambda0 lambda0dot step))
    (multiple-value-bind (dndt-prev dldt-prev d2ndt2-prev d2ldt2-prev)
        (derivs-synchronous-resonance v-deltai lambda0dot n0 lambda0)
      (multiple-value-bind (n l) (integrate-synchronous-resonance
                                  v-deltai lambda0dot
                                  0.0d0 n0 lambda0 dndt-prev dldt-prev d2ndt2-prev d2ldt2-prev
                                  time step)
        (let ((M (- (+ l thetat) (+ O w))))
          (values n M))))))

(defun integrate-from-t50-half-day-resonance
    (t50 v-el0
     v-deldt-Earth-zonal-harmonics
     v-deldt-lunar v-deldt-solar v-Dlmpq
     time v-el)
  "Do Euler-Maclaurin integration for Earth gravity half-day resonance effects starting from epoch
time and going to time TIME minutes after epoch.  V-EL are the elements at time TIME after accounting
for Earth zonal gravity, partial drag and lunar/solar gravity secular effects.  Return (values n M)."
  (declare (type double-float t50 time))
  (let* ((n0 (aref v-el0 0))
         (O0 (aref v-el0 3))
         (w0 (aref v-el0 4))
         (M0 (aref v-el0 5))
         (O (aref v-el 3))
         (dOdt-EZ (aref v-deldt-Earth-zonal-harmonics 3))
         (dwdt-EZ (aref v-deldt-Earth-zonal-harmonics 4))
         (dMdt-EZ (aref v-deldt-Earth-zonal-harmonics 5))
         (dOdt-LS (+ (aref v-deldt-lunar 3) (aref v-deldt-solar 3)))
         ;; (dwdt-LS (+ (aref v-deldt-lunar 4) (aref v-deldt-solar 4)))
         (dMdt-LS (+ (aref v-deldt-lunar 5) (aref v-deldt-solar 5)))
         (dwdt0 dwdt-EZ) ;seems it shouldn't include lunar/solar seculars (???)
         (theta0 (phi-sgp4-at-t50 t50))
         (dthetadt 4.37526908801129966d-3) ;gmst/sidereal rotation rate in radians/minute
         (thetat (+ theta0 (* dthetadt time))) ;should be close enough
         (lambda0 (rdn-mpi-ppi (- (+ M0 (* 2.0d0 O0)) (* 2.0d0 theta0))))
         (lambda0dot (- (+ dMdt-EZ dMdt-LS
                           (* 2.0d0 (+ dOdt-EZ dOdt-LS))) (* 2.0d0 dthetadt)))
         (step 720.0d0))
    (declare (type double-float
                   n0 O0 w0 M0 O dOdt-EZ dwdt-EZ dMdt-EZ dOdt-LS dMdt-LS dwdt0
                   theta0 dthetadt thetat lambda0 lambda0dot step))
    (multiple-value-bind (dndt-prev dldt-prev d2ndt2-prev d2ldt2-prev)
        (derivs-half-day-resonance v-Dlmpq lambda0dot w0 n0 lambda0)
      (multiple-value-bind (n l) (integrate-half-day-resonance
                                  v-Dlmpq w0 dwdt0 lambda0dot
                                  0.0d0 n0 lambda0 dndt-prev dldt-prev d2ndt2-prev d2ldt2-prev
                                  time step)
        (let ((M (- (+ l (* 2.0d0 thetat)) (* 2.0d0 O))))
          (values n M))))))


;;;; SECULAR UPDATE FOR REMAINING ATMOSPHERIC DRAG EFFECTS

(defun update-v-el-for-remaining-drag
    (v-grav v-el0 Bstar v-Ci v-Di drop-terms? time v-el)
  "Return a new elements set, an update from set V-EL accounting for remaining drag effects.
See Hoots2004, App B, section B.4."
  (declare (type double-float Bstar time))
  (let* ((ke (aref v-grav 1))
         (M0 (aref v-el0 5))
         (n (aref v-el 0))
         (e (aref v-el 1))
         (I (aref v-el 2))
         (O (aref v-el 3))
         (w (aref v-el 4))
         (M (aref v-el 5))
         (C1 (aref v-Ci 0))
         (C4 (aref v-Ci 3))
         (C5 (aref v-Ci 4))
         (D2 (aref v-Di 0))
         (D3 (aref v-Di 1))
         (D4 (aref v-Di 2))
         (t2 (* time time))
         (t3 (* time t2))
         (t4 (* t2 t2))
         (e1 (if drop-terms?
                (- e (* Bstar C4 time))
                        ;???? Hoots has e0 but that would seem to mean that
                        ;secular lunar/solar is lost; java program has "e"
                        ;but it seems to be ecco; however, it's applied
                        ;later, so not really clear here
                (- e (+ (* Bstar C4 time)
                        (* Bstar C5 (- (sin M) (sin M0)))))))
         (e2 (if (< e1 1.0d-6) 1.0d-6 e1))
         (a (let ((bkt (if drop-terms?
                           (- 1.0d0 (* C1 time))
                           (- 1.0d0 (+ (* C1 time) (* D2 t2)
                                       (* D3 t3) (* D4 t4))))))
              (* (expt (/ ke n) (/ 2.0d0 3.0d0)) (* bkt bkt))))
         (n (/ ke (expt a 1.5d0))))
    (declare (type double-float ke M0 n e e1 e2 I O w M C1 C4 C5 D2 D3 D4 t2 t3 t4 e e a n))
    (make-v-el n e2 I O w M)))


;;;; IL CALCULATION

(defun calculate-IL
    (v-el0 v-Ci v-Di drop-terms? time v-el)
  "Calculate IL value of Hoots, App B, section B.4."
  (declare (type double-float time))
  (let* ((n0 (aref v-el0 0))
         (O (aref v-el 3))
         (w (aref v-el 4))
         (M (aref v-el 5))
         (C1 (aref v-Ci 0))
         (D2 (aref v-Di 0))
         (D3 (aref v-Di 1))
         (D4 (aref v-Di 2))
         (t2 (* time time))
         (t3 (* time t2))
         (t4 (* t2 t2))
         (t5 (* time t4))
         (C1to2 (* C1 C1))
         (C1to3 (* C1 C1to2))
         (C1to4 (* C1to2 C1to2))
         (D2to2 (* D2 D2)))
    (declare (type double-float n0 O w M C1 D2 D3 D4 t2 t3 t4 t5 C1to2 C1to3 C1to4 D2to2))
    ;; could work out the coefficients as part of initialization
    (if drop-terms?
        (+ M w O (* n0 (* 1.5d0 C1 t2)))
        (+ M w O
           (* n0 (+ (* 1.5d0 C1 t2) (* (+ D2 (* 2.0d0 C1to2)) t3)
                    (* 0.25d0 (+ (* 3.0d0 D3) (* 12.0d0 C1 D2)
                                   (* 10.0d0 C1to3)) t4)
                    (* 0.2d0 (+ (* 3.0d0 D4) (* 12.0d0 C1 D3)
                                  (* 6.0d0 D2to2) (* 30.0d0 C1to2 D2)
                                  (* 15.0d0 C1to4)) t5)))))))


;;;; UPDATE FOR LONG PERIOD PERIODIC EFFECTS OF LUNAR AND SOLAR GRAVITY

(defun v-dellp-for-long-period-perturber
    (v-el0 v-Xi v-Zi v-Zij ex Mx Cx)
  "Return a 5-vector [de dI dM dwpcIdO sIdO] for long period lunar/solar perturbations.  Here de, dI
and dM are deltas for eccentricity inclination and mean anomaly, respectively.  dwpcIdO
is (\\delta\\omega + \\cos I \\delta\\Omega).  sIdO is (\\sin I \\delta\\Omega).  Input V-EL0 are the
satellite epoch elements.  Inputs 'v-Xi', 'v-Zi' and 'v-Zij' are from an sgp4data record for the
perturbing body.  Inputs 'ex', 'Mx' and 'Cx' are the perturbing body's orbital eccentricity, mean
anomaly and perturbation coefficient respectively.  See Hoots2004, App B, section B.5 and Hoots2004,
App A section E."
  (declare (type double-float ex Mx Cx))
  (let* ((n0 (aref v-el0 0))
         (e0 (aref v-el0 1))
         (X1 (aref v-Xi 0))
         (X2 (aref v-Xi 1))
         (X3 (aref v-Xi 2))
         (X4 (aref v-Xi 3))
         (Z1 (aref v-Zi 0))
         (Z2 (aref v-Zi 1))
         (Z3 (aref v-Zi 2))
         (Z11 (aref v-Zij 0))
         (Z12 (aref v-Zij 1))
         (Z13 (aref v-Zij 2))
         (Z21 (aref v-Zij 3))
         (Z22 (aref v-Zij 4))
         (Z23 (aref v-Zij 5))
         (Z31 (aref v-Zij 6))
         (Z32 (aref v-Zij 7))
         (Z33 (aref v-Zij 8))
         ;; (Mx (rdn-0-2pi Mx)) ;no effect on 23333
         (fx (+ Mx (* 2.0d0 ex (sin Mx))))
         (cosfx (cos fx))
         (sinfx (sin fx))
         (F2 (- (* 0.5d0 (* sinfx sinfx)) 0.25d0))
         (F3 (- (* 0.5d0 sinfx cosfx)))
         (eta0 (sqrt (- 1.0d0 (* e0 e0)))) ;???? seems not defined in Hoots2004
         ;; "lpcoeff-xx-yy" is "long period coefficient" in the formula for
         ;; xx, multiplying yy
         ;; could work out coefficients during initialization
         (lpcoeff-de-part (/ (* -30.0d0 eta0 Cx e0) n0))
         (lpcoeff-de-F2 (* lpcoeff-de-part (+ (* X2 X3) (* X1 X4))))
         (lpcoeff-de-F3 (* lpcoeff-de-part (- (* X2 X4) (* X1 X3))))
         (de (+ (* lpcoeff-de-F2 F2) (* lpcoeff-de-F3 F3)))
         (lpcoeff-dI-part (- (/ Cx (* n0 eta0))))
         (lpcoeff-dI-F2 (* lpcoeff-dI-part Z12))
         (lpcoeff-dI-F3 (* lpcoeff-dI-part (- Z13 Z11)))
         (dI (+ (* lpcoeff-dI-F2 F2) (* lpcoeff-dI-F3 F3)))
         (lpcoeff-dM-part (/ (* -2.0d0 Cx) n0))
         (lpcoeff-dM-F2 (* lpcoeff-dM-part Z2))
         (lpcoeff-dM-F3 (* lpcoeff-dM-part (- Z3 Z1)))
         (lpcoeff-dM-sinfx (* -3.0d0 lpcoeff-dM-part ex
                              (+ 7.0d0 (* 3.0d0 (* e0 e0)))))
                        ;???? Hoots has e_x not e_X ????
         (dM (+ (* lpcoeff-dM-F2 F2) (* lpcoeff-dM-F3 F3)
                (* lpcoeff-dM-sinfx sinfx)))
         (lpcoeff-dwpcIdO-part (/ (* 2.0d0 eta0 Cx) n0))
         (lpcoeff-dwpcIdO-F2 (* lpcoeff-dwpcIdO-part Z32))
         (lpcoeff-dwpcIdO-F3 (* lpcoeff-dwpcIdO-part (- Z33 Z31)))
         (lpcoeff-dwpcIdO-sinfx (* lpcoeff-dwpcIdO-part (* -9.0d0 ex)))
         (dwpcIdO (+ (* lpcoeff-dwpcIdO-F2 F2) (* lpcoeff-dwpcIdO-F3 F3)
                     (* lpcoeff-dwpcIdO-sinfx sinfx)))
                        ;\delta \omega_x + \cos I_x \delta \Omega_x
         (lpcoeff-sIdO-part (/ Cx (* n0 eta0)))
         (lpcoeff-sIdO-F2 (* lpcoeff-sIdO-part Z22))
         (lpcoeff-sIdO-F3 (* lpcoeff-sIdO-part (- Z23 Z21)))
         (sIdO (+ (* lpcoeff-sIdO-F2 F2) (* lpcoeff-sIdO-F3 F3)))
                        ;\sin I_x \delta \Omega_x
         )
    (declare (type double-float
                   n0 e0 X1 X2 X3 X4 Z1 Z2 Z3 Z11 Z12 Z13 Z21 Z22 Z23 Z31 Z32 Z33
                   fx cosfx sinfx F2 F3 eta0 lpcoeff-de-part lpcoeff-de-F2 lpcoeff-de-F3
                   de lpcoeff-dI-part lpcoeff-dI-F2 lpcoeff-dI-F3 dI lpcoeff-dM-part lpcoeff-dM-F2
                   lpcoeff-dM-F3 lpcoeff-dM-sinfx dM lpcoeff-dwpcIdO-part
                   lpcoeff-dwpcIdO-F2 lpcoeff-dwpcIdO-F3 lpcoeff-dwpcIdO-sinfx dwpcIdO
                   lpcoeff-sIdO-part lpcoeff-sIdO-F2 lpcoeff-sIdO-F3 sIdO))
    ;; (format t "Mx, ex, fx sinfx cosfx = ~a ~a ~a ~a ~a~%" Mx ex fx sinfx cosfx)
    ;; (format t "F2, F3 = ~a ~a~%" F2 F3)
    ;; (format t "de, dI, dM = ~a ~a ~a~%" de dI dM)
    ;; (format t "dwpcIdO, sIdO = ~a ~a~%" dwpcIdO sIdO)
    ;; lunar ones match printout of SEL, SIL, SLL, SGHL, SHL,
    ;; solar of SES, SIS, SLS, SGHS, SHS,
    ;; but not too accurately, say to 1e-8 or so
    ;; lunar F2, F3 seem worse than solar
    (vector de dI dM dwpcIdO sIdO)))

(defun v-dellp-for-long-period-solar (t50 v-el0 v-Xi v-Zi v-Zij time)
  "Return a 5-vector [de dI dM dwpcIdO sIdO] for long period solar perturbations (see docstring of
v-dellp-for-long-period-perturber for a description of these quantities).  Inputs T50 and V-EL0 are
the satellite epoch time and epoch elements.  Inputs 'v-Xi', 'v-Zi' and 'v-Zij' are from an sgp4data
record for solar perturbations.  Inputs TIME and V-EL are minutes after epoch and orbital elements of
the satellite at time TIME respectively."
  (declare (type double-float t50 time))
  (let* ((ns 1.19459d-5) ;solar mean motion
         (es 0.01675d0)    ;solar eccentricity
         ;;Is  0.4091767351668026d0 ;radians = 23.4441 deg
         (Cs 2.9864797d-6) ;C_s solar perurbation coefficient, radians/min
         (day (+ t50 18261.5d0))
         (Ms0 (rdn-mpi-ppi (+ 6.2565837d0 (* 0.017201977d0 day)))) ;at epoch
         (Ms (+ Ms0 (* ns time))))
    (declare (type double-float ns es Cs day Ms0 Ms))
    ;; (format t "solar time = ~a~%" time)
    (v-dellp-for-long-period-perturber v-el0 v-Xi v-Zi v-Zij es Ms Cs)))

(defun v-dellp-for-long-period-lunar
    (t50 v-el0 v-Xi v-Zi v-Zij time)
  "Return a 5-vector [de dI dM dwpcIdO sIdO] for long period lunar perturbations (see docstring of
v-dellp-for-long-period-perturber for a description of these quantities).  Inputs T50 and V-EL0 are
the satellite epoch time and epoch elements.  Inputs 'v-Xi', 'v-Zi' and 'v-Zij' are from an sgp4data
record for lunar perturbations.  Inputs TIME and V-EL are minutes after epoch and orbital elements of
the satellite at time TIME respectively."
  (declare (type double-float t50 time))
  (let* (;; nl 1.583521770d-4 ;lunar mean motion from Hoots2004
         (nl 1.5835218d-4) ;lunar mean motion from existing program
         (el 0.05490d0)      ;lunar eccentricity
         (Cl 4.796806521d-7) ;C_m lunar perurbation coefficient, radians/min
         (day      (+ t50 18261.5d0)) ;lunar ephemeris uses Jan 0.5, 1900.0
         ;; (Ome      (rdn-mpi-ppi (- 4.5236020 (* 9.2422029d-4 day))))
                        ;\Omega_m_epsilon (?)
         ;; (cosOme   (cos Ome)) ;cos(\Omega_m_\epsilon)
         ;; (sinOme   (sin Ome)) ;sin(\Omega_m_\epsilon)
         ;; (cosIm    (- 0.91375164d0 (* 0.03568096d0 cosOme)))
                        ;cos(I_m); here 0.913... = cos(\epsilon)
                        ;cos(I_m_\epsilon); 0.035... = sin(\epsilon)
                        ;sin(I_m_\epsilon)
         (gamma    (+ 5.8351514d0 (* 0.0019443680d0 day)))
                        ;lunar longitude of perigee referred to the ecliptic
         ;;(Ml0 (- (+ 4.7199672d0 (* 0.22997150d0 day)) gamma))
         (Ml0 (rdn-mpi-ppi (- (+ 4.7199672d0 (* 0.22997150d0 day)) gamma)))
         (Ml (+ Ml0 (* nl time))))
    (declare (type double-float nl el Cl day gamma Ml0 Ml))
    ;; (format t "lunar time = ~a~%" time)
    (v-dellp-for-long-period-perturber v-el0 v-Xi v-Zi v-Zij el Ml Cl)))

(defun update-v-el-for-long-period-perturber (v-el v-dellp)
  "Return a 6-vector element set which is input set V-EL updated for the long period perturbations
defined by the element deltas 'v-dellp'.  Lyddane action is all within this function."
  (let* ((n (aref v-el 0))
         (e (aref v-el 1))
         (I (aref v-el 2))
         (O (aref v-el 3))
         (w (aref v-el 4))
         (M (aref v-el 5))
         (O (rdn-0-2pi O)) ;without this, test case 23599 (which has O close to 0) is wrong by 960 m
         (de (aref v-dellp 0))
         (dI (aref v-dellp 1))
         (dM (aref v-dellp 2))
         (dwpcIdO (aref v-dellp 3))
         (sIdO (aref v-dellp 4))
         (e (+ e de))
         (I (+ I dI))
         (sinI (sin I))
         (cosI (cos I))
         (lyddane? (<= I 0.2)) ;Vallado's "GSFC choice"
         (dO (if lyddane? 0.0d0 (/ sIdO sinI))) ;not used if lyddane? is true
         (O-old O)
         (O (if lyddane?
                (let* ((sinO (sin O))
                       (cosO (cos O))
                       (alpha (+ (* sinI sinO) (* cosO sIdO)
                                 (* cosI sinO dI)))
                       (beta (+ (- (* sinI cosO) (* sinO sIdO))
                                (* cosI cosO dI)))
                       (O-1 (atan alpha beta))
                       (O-2 (if (> (abs (- O O-1)) pi)
                                (if (< O-1 O) (+ O-1 (* 2.0d0 pi)) (- O-1 (* 2.0d0 pi)))
                                O-1)))
                  O-2)
                (+ O dO) ;not what it says to do in Hoots as far as I can
                        ;see, but matches better with the canonical program
                ))
         (w (if lyddane?
                ;; (+ w (* O-old (cos I)) dwpcIdO
                ;; (- (* dI O-old (sin I))) (- (* O (cos I))))
                ;; ???? this is my best stab at what the existing program does
                ;; it simplifies further to:
                (- (+ w (- dwpcIdO (* cosI (- O O-old)))) (* dI sinI O-old))
                (+ w (- dwpcIdO (* cosI dO)))
                ;; ???? not what it says to do in Hoots as far as I
                ;; can see, but matches better with the program
                ))
         (M (+ M dM)))
    (declare (type double-float n e I O w M O de dI dM dwpcIdO sIdO e I sinI cosI dO O-old O w M))
    (make-v-el-fixing-negative-I n e I O w M)))

(defun update-for-long-period-lunar-solar
    (t50 v-el0
     v-Xi-lunar v-Zi-lunar v-Zij-lunar
     v-Xi-solar v-Zi-solar v-Zij-solar
     time v-el)
  "Return a 6-vector element set which is input set V-EL updated for the long period lunar and solar
perturbations.  Other inputs are from an sgp4data record for the lunar and solar perturbations."
  (declare (type double-float t50 time))
  (let* ((v-dellp-solar (v-dellp-for-long-period-solar
                         t50 v-el0 v-Xi-solar v-Zi-solar v-Zij-solar time))
         (v-dellp-lunar (v-dellp-for-long-period-lunar
                         t50 v-el0 v-Xi-lunar v-Zi-lunar v-Zij-lunar time))
         (v-dellp (make-array 5 :element-type 'double-float)))
    (loop for i below 5
       do (setf (aref v-dellp i) (+ (aref v-dellp-solar i) (aref v-dellp-lunar i))))
    (update-v-el-for-long-period-perturber v-el v-dellp)))


;;;;; UPDATE FOR LONG PERIOD PERIODIC EFFECTS OF EARTH GRAVITY

(defun update-for-long-period-Earth-gravity (v-grav v-el IL)
  "Return a (values axN ILT ayN), values needed to update to account for long period periodic effects
of Earth gravity.  See Hoots2004, App B, section B.6."
  (declare (type double-float IL))
  (let* ((ke (aref v-grav 1))
         (J2 (aref v-grav 2))
         (J3 (aref v-grav 3))
         (k2 (* 0.5d0 J2)) ;{1 \over 2} J_2 a_E^2 but a_E in Earth radii is 1
         (A30 (- J3)) ;-J_3 a_E^3 but a_E in Earth radii is 1
         (n (aref v-el 0))
         (e (aref v-el 1))
         (I (aref v-el 2))
         (w (aref v-el 4))
         (a (expt (/ ke n) (/ 2.0d0 3.0d0)))
         (beta2 (- 1.0d0 (* e e))) ;\beta^2
         (cosI (cos I))
         (sinI (sin I))
         (cosw (cos w))
         (sinw (sin w))
         (axN (* e cosw))
         (ILL (* (/ (* A30 sinI) (* 8.0d0 k2 a beta2))
                 axN    ;==(* e cosw)
                 (/ (+ 3.0d0 (* 5.0d0 cosI)) (max (+ 1.0d0 cosI) 1.5d-12))))
                        ;fix for divide by zero with I = 180 deg
         (ayNL (/ (* A30 sinI) (* 4.0d0 k2 a beta2)))
         (ILT (+ IL ILL))
         (ayN (+ (* e sinw) ayNL)))
    (declare (type double-float
                   ke J2 J3 k2 A30 n e I w a beta2 cosI sinI cosw sinw axN ILL ayNL ILT ayN))
    (values axN ILT ayN)))


;;;; UPDATE FOR SHORT PERIOD PERIODIC EFFECTS OF EARTH GRAVITY

(defun update-for-short-period-Earth-gravity (v-grav v-el axN ILT ayN)
  "Return (values rk uk Ok Ik drkdt rdfkdt), values needed to update to account for short period
periodic effects of Earth gravity.  See Hoots2004, App B, section B.7."
  (declare (type double-float axN ILT ayN))
  (let* ((ke (aref v-grav 1))
         (J2 (aref v-grav 2))
         (k2 (* 0.5d0 J2)) ;{1 \over 2} J_2 a_E^2 but a_E in Earth radii is 1
         (n (aref v-el 0))
         (I (aref v-el 2))
         (O (aref v-el 3))
         ;; reduction to 0-2pi or mpi-ppi not in Hoots, but is in the program
         (ILT (rdn-mpi-ppi ILT))
         (O (rdn-mpi-ppi O))
         (U (- ILT O))
         (Epw (loop     ;solve Kepler's equation; "Epw" = E+\omega
                 for iter of-type fixnum from 0
                 with Epw of-type double-float = (coerce U 'double-float)
                 do (let* ((cosEpw (cos Epw))
                           (sinEpw (sin Epw))
                           (DEpw-1 (/ (- (+ U (* axN sinEpw)) (+ (* ayN cosEpw) Epw))
                                      (- 1.0d0 (+ (* ayN sinEpw) (* axN cosEpw)))))
                           (DEpw (if (> (abs DEpw-1) 1.0d0) ;in the Fortran program, but not in Hoots
                                     (float-sign DEpw-1)
                                     DEpw-1)))
                      (declare (type double-float cosEpw sinEpw DEpw-1 DEpw))
                      (cond
                        ((< (abs DEpw) 1.0d-12) (return (+ Epw DEpw)))
                        ((>= iter 20)
                         (error 'Keplers-equation-convergence-error :text
                                (format nil "Kepler's eqn iteration limit exceeded, v-el = ~%~a~%"
                                        v-el))
                         ;; if continued (proper handling not done yet):
                         (return (+ Epw DEpw)))
                        (t (setf Epw (+ Epw DEpw)))))))
         (cosEpw (cos Epw))
         (sinEpw (sin Epw))
         (ecosE (+ (* axN cosEpw) (* ayN sinEpw)))
         (esinE (- (* axN sinEpw) (* ayN cosEpw)))
         (e2 (+ (* axN axN) (* ayN ayN)))
         (oneme2 (- 1.0d0 e2))
         (a (expt (/ ke n) (/ 2.0d0 3.0d0)))
         (pL (* a oneme2))
         (r (* a (- 1.0d0 ecosE)))
         (drdt (* ke (/ (sqrt a) r) esinE))
         (rdfdt (* ke (/ (sqrt pL) r)))
         (aoverr (/ a r))
         (onepsqrtoneme2 (+ 1.0d0 (sqrt oneme2)))
         (cosu (* aoverr (+ (- cosEpw axN) (/ (* ayN esinE) onepsqrtoneme2))))
         (sinu (* aoverr (- sinEpw ayN (/ (* axN esinE) onepsqrtoneme2))))
         (u (atan sinu cosu))
         (cosI (cos I)) ;;???? note that STR3 uses cos I0 here
         (sinI (sin I))
         (cosI2 (* cosI cosI))
         (cos2u (cos (* 2.0d0 u)))
         (sin2u (sin (* 2.0d0 u)))
         (pL2 (* pL pL))
         (DO-DI-common (/ (* 3.0d0 k2 cosI) (* 2.0d0 pL2)))
         (Dr (* (/ k2 (* 2.0d0 pL)) (- 1.0d0 cosI2) cos2u))
         (Du (* (/ k2 (* -4.0d0 pL2)) (- (* 7.0d0 cosI2) 1.0d0) sin2u))
         (DO (* DO-DI-common sin2u))
         (DI (* DO-DI-common sinI cos2u))
         (Ddrdt (* (- (/ (* k2 n) pL)) (- 1.0d0 cosI2) sin2u))
         (Drdfdt (* (/ (* k2 n) pL) (- (* (- 1.0d0 cosI2) cos2u)
                                       (* 1.5d0 (- 1.0d0 (* 3.0d0 cosI2))))))
         (rk (+ (* r (- 1.0d0 (* 1.5d0 k2 (/ (sqrt oneme2) pL2)
                                 (- (* 3.0d0 cosI2) 1.0d0)))) Dr))
         (uk (+ u Du))
         (Ok (+ O DO))
         (Ik (+ I DI))
         (drkdt (+ drdt Ddrdt))
         (rdfkdt (+ rdfdt Drdfdt)))
    (declare (type double-float
                   ke J2 k2 n I O ILT U Epw cosEpw sinEpw ecosE esinE e2 oneme2 a pL r
                   drdt rdfdt aoverr onepsqrtoneme2 cosu sinu u cosI sinI cosI2 cos2u sin2u pL2 Dr Du
                   DO DI Ddrdt Drdfdt rk uk Ok Ik drkdt rdfkdt))
    (values rk uk Ok Ik drkdt rdfkdt)))


;;;;; MAKE FINAL V-STATE IN CANONICAL UNITS

(defun make-final-v-state-canonical-units (rk uk Ok Ik drkdt rdfkdt)
  "Return a 6-vector containing [x y z vx vy vz], the position and velocity state in SGP4 canonical
units (Earth radii and Earth radii per minute), given the arguments which are data as generated by
update-for-short-period-Earth-gravity."
  (declare (type double-float rk uk Ok Ik drkdt rdfkdt))
  (let* (;; orientation vectors
         (cosOk (cos Ok))
         (sinOk (sin Ok))
         (cosIk (cos Ik))
         (sinIk (sin Ik))
         (cosuk (cos uk))
         (sinuk (sin uk))
         (M0 (- (* sinOk cosIk)))
         (M1  (* cosOk cosIk))
         (M2  sinIk)
         (N0 cosOk)
         (N1 sinOk)
         (N2 0.0d0)
         (U0 (+ (* M0 sinuk) (* N0 cosuk)))
         (U1 (+ (* M1 sinuk) (* N1 cosuk)))
         (U2 (+ (* M2 sinuk) (* N2 cosuk)))
         (V0 (- (* M0 cosuk) (* N0 sinuk)))
         (V1 (- (* M1 cosuk) (* N1 sinuk)))
         (V2 (- (* M2 cosuk) (* N2 sinuk)))
         (r0 (* rk U0))
         (r1 (* rk U1))
         (r2 (* rk U2))
         (v0 (+ (* drkdt U0) (* rdfkdt V0)))
         (v1 (+ (* drkdt U1) (* rdfkdt V1)))
         (v2 (+ (* drkdt U2) (* rdfkdt V2)))
         (v-state (make-array 6 :element-type 'double-float)))
    (declare (type double-float cosOk sinOk cosIk sinIk cosuk sinuk
                   M0 M1 M2 N0 N1 N2 U0 U1 U2 V0 V1 V2 r0 r1 r2 v0 v1 v2))
    (when (< (+ (* r0 r0) (* r1 r1) (* r2 r2)) 1.0d0)
      (error 'altitude-lt-0-error :text "altitude < 0"))
    (setf (aref v-state 0) r0)
    (setf (aref v-state 1) r1)
    (setf (aref v-state 2) r2)
    (setf (aref v-state 3) v0)
    (setf (aref v-state 4) v1)
    (setf (aref v-state 5) v2)
    v-state))


;;;; USER LEVEL PROPAGATION FUNCTION

(defun propagate-sgp4 (sgp4data time)
  "Return state vector, [x y z vx vy vz] in canonical units (Earth radii and minutes), at TIME minutes
after epoch.  'sgp4data' is an SGP4 initialization record as generated by function init-sgp4."
  (declare (type sgp4data sgp4data) (type double-float time))
  (let* (;; ! clos class would be better for sgp4data just to be able to use with-slots !
         (v-grav (sgp4data-v-grav sgp4data))
         (t50 (sgp4data-t50 sgp4data))
         (v-el0 (sgp4data-v-el0 sgp4data))
         (Bstar (sgp4data-Bstar sgp4data))
         (q0 (sgp4data-q0 sgp4data))
         (s (sgp4data-s sgp4data))
         (v-Ci (sgp4data-v-Ci sgp4data))
         (v-Di (sgp4data-v-Di sgp4data))
         (v-deldt-Earth-zonal-harmonics (sgp4data-v-deldt-Earth-zonal-harmonics sgp4data))
         (deep-space? (sgp4data-deep-space? sgp4data))
         (resonance (sgp4data-resonance sgp4data))
         (drop-terms? (sgp4data-drop-terms? sgp4data))
         (v-Xi-solar (sgp4data-v-Xi-solar sgp4data))
         (v-Zi-solar (sgp4data-v-Zi-solar sgp4data))
         (v-Zij-solar (sgp4data-v-Zij-solar sgp4data))
         (v-deldt-solar (sgp4data-v-deldt-solar sgp4data))
         (v-Xi-lunar (sgp4data-v-Xi-lunar sgp4data))
         (v-Zi-lunar (sgp4data-v-Zi-lunar sgp4data))
         (v-Zij-lunar (sgp4data-v-Zij-lunar sgp4data))
         (v-deldt-lunar (sgp4data-v-deldt-lunar sgp4data))
         (v-deltai (sgp4data-v-deltai sgp4data))
         (v-Dlmpq (sgp4data-v-Dlmpq sgp4data))
         (v-el1 (update-v-el-for-Earth-zonal-gravity-and-partial-drag
                 v-el0 v-grav Bstar q0 s v-Ci
                 v-deldt-Earth-zonal-harmonics drop-terms?
                 time))
         (v-el2 (if deep-space?
                    (update-v-el-for-secular-lunar-solar-gravity
                     v-el1 v-deldt-lunar v-deldt-solar time)
                    v-el1))
         (v-el3 (case resonance
                        ;integrate from epoch time (could add state)
                  (:synchronous
                   (multiple-value-bind (n M)
                       (integrate-from-t50-synchronous-resonance
                        t50 v-el0 v-deldt-Earth-zonal-harmonics
                        v-deldt-lunar v-deldt-solar v-deltai
                        time v-el2)
                     (copy-v-el v-el2  0 n  5 M)))
                  (:half-day
                   (multiple-value-bind (n M)
                       (integrate-from-t50-half-day-resonance
                        t50 v-el0 v-deldt-Earth-zonal-harmonics
                        v-deldt-lunar v-deldt-solar v-Dlmpq
                        time v-el2)
                     (copy-v-el v-el2  0 n  5 M)))
                  ;; else, no resonance
                  (otherwise v-el2)))
         (v-el4 (update-v-el-for-remaining-drag
                 v-grav v-el0 Bstar v-Ci v-Di drop-terms? time v-el3))
         (v-el5 (if deep-space?
                    (update-for-long-period-lunar-solar
                     t50 v-el0
                     v-Xi-lunar v-Zi-lunar v-Zij-lunar
                     v-Xi-solar v-Zi-solar v-Zij-solar
                     time v-el4)
                    v-el4))
         (IL (calculate-IL v-el0 v-Ci v-Di drop-terms? time v-el5)))
    (multiple-value-bind (axN ILT ayN) (update-for-long-period-Earth-gravity v-grav v-el5 IL)
      (multiple-value-bind (rk uk Ok Ik drkdt rdfkdt)
          (update-for-short-period-Earth-gravity v-grav v-el5 axN ILT ayN)
        (make-final-v-state-canonical-units rk uk Ok Ik drkdt rdfkdt)))))


;;;; STATE VECTOR UNIT CONVERSIONS

;; NOTE --
;; MAKE-V-GRAV-WGS-72 and 72-LOW-PRECISION use the value for Earth radius given here;
;; MAKE-V-GRAV-WGS-84 uses 6378.137

(defparameter *km-per-unit* 6378.135d0)
(defparameter *km/s-per-unit* (/ 6378.135d0 60.0d0))

(defparameter *m-per-unit* 6378135.0d0)
(defparameter *m/s-per-unit* (/ 6378135.0d0 60.0d0))

(defmacro with-v-grav-units (v-grav &rest body)
  "Execute BODY with unit converters appropriate for the Earth radius specified in V-GRAV.  V-GRAV is
the vector of gravity parameters in the form returned by one of the MAKE-V-GRAV-xxx functions.  The
unit converters affected are the functions ->v-state-km-s and ->v-state-m-s which convert canonical
units to km, km/s and m, m/s respectively."
  (let ((cur-km-per-unit (gensym))
        (cur-km/s-per-unit (gensym))
        (cur-m-per-unit (gensym))
        (cur-m/s-per-unit (gensym))
        (val (gensym)))
    `(let ((,cur-km-per-unit *km-per-unit*)
           (,cur-km/s-per-unit *km/s-per-unit*)
           (,cur-m-per-unit *m-per-unit*)
           (,cur-m/s-per-unit *m/s-per-unit*))
       (setf *km-per-unit* (aref ,v-grav 0))
       (setf *km/s-per-unit* (/ (aref ,v-grav 0) 60.0d0))
       (setf *m-per-unit* (* 1000.0d0 (aref ,v-grav 0)))
       (setf *m/s-per-unit* (/ (* 1000.0d0 (aref ,v-grav 0)) 60.0d0))
       (let ((,val ,@body))
         (setf *km-per-unit* ,cur-km-per-unit)
         (setf *km/s-per-unit* ,cur-km/s-per-unit)
         (setf *m-per-unit* ,cur-m-per-unit)
         (setf *m/s-per-unit* ,cur-m/s-per-unit)
         ,val))))

(defun ->v-state-km-s (v-state)
  "Convert input state vector V-STATE, a 6-vector containing [x y z vx vy vz] in canonical units,
to a state 6-vector in km, km/s units."
  (declare (type (array double-float (6)) v-state))
  (let* ((v-state-units (make-array 6 :element-type 'double-float)))
    (setf (aref v-state-units 0) (* *km-per-unit* (aref v-state 0)))
    (setf (aref v-state-units 1) (* *km-per-unit* (aref v-state 1)))
    (setf (aref v-state-units 2) (* *km-per-unit* (aref v-state 2)))
    (setf (aref v-state-units 3) (* *km/s-per-unit* (aref v-state 3)))
    (setf (aref v-state-units 4) (* *km/s-per-unit* (aref v-state 4)))
    (setf (aref v-state-units 5) (* *km/s-per-unit* (aref v-state 5)))
    v-state-units))

(defun ->v-location-km (v-state)
  "Return a 3-vector containing the satellite location in km corresponding to input state vector
V-STATE.  V-STATE is a 6-vector containing [x y z vx vy vz] in canonical units, output from
PROPAGATE-SGP4,"
  (declare (type (array double-float (6)) v-state))
  (let* ((v-location-km (make-array 3 :element-type 'double-float)))
    (setf (aref v-location-km 0) (* *km-per-unit* (aref v-state 0)))
    (setf (aref v-location-km 1) (* *km-per-unit* (aref v-state 1)))
    (setf (aref v-location-km 2) (* *km-per-unit* (aref v-state 2)))
    v-location-km))

(defun ->v-velocity-km/s (v-state)
  "Return a 3-vector containing the satellite velocity in km/s corresponding to input state vector
V-STATE.  V-STATE is a 6-vector containing [x y z vx vy vz] in canonical units, output from
PROPAGATE-SGP4,"
  (declare (type (array double-float (6)) v-state))
  (let* ((v-velocity-km/s (make-array 3 :element-type 'double-float)))
    (setf (aref v-velocity-km/s 0) (* *km/s-per-unit* (aref v-state 3)))
    (setf (aref v-velocity-km/s 1) (* *km/s-per-unit* (aref v-state 4)))
    (setf (aref v-velocity-km/s 2) (* *km/s-per-unit* (aref v-state 5)))
    v-velocity-km/s))

(defun ->v-state-m-s (v-state)
  "Convert input state vector V-STATE, a 6-vector containing [x y z vx vy vz] in canonical units,
to a state 6-vector in m, m/s units."
  (declare (type (array double-float (6)) v-state))
  (let* ((v-state-units (make-array 6 :element-type 'double-float)))
    (setf (aref v-state-units 0) (* *m-per-unit* (aref v-state 0)))
    (setf (aref v-state-units 1) (* *m-per-unit* (aref v-state 1)))
    (setf (aref v-state-units 2) (* *m-per-unit* (aref v-state 2)))
    (setf (aref v-state-units 3) (* *m/s-per-unit* (aref v-state 3)))
    (setf (aref v-state-units 4) (* *m/s-per-unit* (aref v-state 4)))
    (setf (aref v-state-units 5) (* *m/s-per-unit* (aref v-state 5)))
    v-state-units))

(defun ->v-location-m (v-state)
  "Return a 3-vector containing the satellite location in meters corresponding to input state
vector V-STATE.  V-STATE is a 6-vector containing [x y z vx vy vz] in canonical units, output
from PROPAGATE-SGP4,"
  (declare (type (array double-float (6)) v-state))
  (let* ((v-location-m (make-array 3 :element-type 'double-float)))
    (setf (aref v-location-m 0) (* *m-per-unit* (aref v-state 0)))
    (setf (aref v-location-m 1) (* *m-per-unit* (aref v-state 1)))
    (setf (aref v-location-m 2) (* *m-per-unit* (aref v-state 2)))
    v-location-m))

(defun ->v-velocity-m/s (v-state)
  "Return a 3-vector containing the satellite velocity in m/s corresponding to input state vector
V-STATE.  V-STATE is a 6-vector containing [x y z vx vy vz] in canonical units, output from
PROPAGATE-SGP4,"
  (declare (type (array double-float (6)) v-state))
  (let* ((v-velocity-m/s (make-array 3 :element-type 'double-float)))
    (setf (aref v-velocity-m/s 0) (* *m/s-per-unit* (aref v-state 3)))
    (setf (aref v-velocity-m/s 1) (* *m/s-per-unit* (aref v-state 4)))
    (setf (aref v-velocity-m/s 2) (* *m/s-per-unit* (aref v-state 5)))
    v-velocity-m/s))
