;; -*- fill-column:102; comment-column:24; -*-

(defpackage :cl-sgp4
  (:use :common-lisp)

  (:documentation
   "Common Lisp implementation of the SGP4 satellite orbit propagation model")

  (:export

   ;; ANGLE CONVERSIONS
   #:deg->rdn
   #:rdn->deg
   #:ndms->rdn
   #:rnd->dms
   #:rdn->ndms
   #:rdn-0-2pi
   #:rdn-mpi-ppi

   ;; SPHERICAL GEOMETRY
   #:v->phi-theta
   #:phi-theta->normalized-v
   #:normalized-v
   #:rdn-between-normalized-v
   #:rdn-between-v
   #:rdn-between-phi-theta

   ;; TIME CONVERSIONS
   #:leap-year-p
   #:nday-in-iyr-imoy
   #:ymd->iJD
   #:iJD->ymd
   #:hms->day
   #:day->hms
   #:day->nhms
   #:ymdhmsZ->JD
   #:iJD->ymd
   #:hms->day
   #:day->hms
   #:day->nhms
   #:ymdhmsZ->JD
   #:JD->ymdhmsZ
   #:t50->JD
   #:JD->t50
   #:ymdhmsZ->t50
   #:t50->ymdhmsZ
   #:universal-time->JD
   #:JD->universal-time
   #:universal-time->t50
   #:t50->universal-time

   ;; SIDEREAL TIME
   #:phi-IAU-at-JD
   #:phi-sgp4-at-t50
   #:phi-IAU-at-t50
   #:phi-IAU-at-universal-time
   #:hms->ra-rdn

   ;; EARTH GEOMETRY
   #:xy-z->latdetic-altitude
   #:txyz->lla


   ;; TLEDATA STRUCTURE
   #:tledata
   #:tledata-num
   #:tledata-classification
   #:tledata-designator-launch-year
   #:tledata-designator-num-of-year
   #:tledata-designator-piece-of-launch
   #:tledata-iyr-norad-epoch ;2-digit year
   #:tledata-doy-norad-epoch
   #:tledata-rev-per-day2-dmm-dt-over-2
   #:tledata-rev-per-day3-dmm-dt2-over-6
   #:tledata-Bstar
   #:tledata-ephermeris-type
   #:tledata-element-num
   #:tledata-line1-checksum
   #:tledata-incl-deg
   #:tledata-ra-of-asc-node-deg
   #:tledata-ecc
   #:tledata-arg-perigee-deg
   #:tledata-mean-anomaly-deg
   #:tledata-mm-rev-per-day
   #:tledata-rev-at-epoch
   #:tledata-line2-checksum

   ;; MAKING TLEDATA STRUCTURES FROM TLE ELEMENT SETS
   #:make-tledata-from-lines
   #:read-tledata-file

   ;; EXTRACTING INFORMATION NEEDED BY SGP4 FROM TLEDATA
   #:tledata->t50-v-el0-Bstar


   ;; SGP4DATA STRUCTURE
   #:sgp4data
   #:sgp4data-v-grav
   #:sgp4data-t50
   #:sgp4data-v-el0
   #:sgp4data-Bstar
   #:sgp4data-q0
   #:sgp4data-s
   #:sgp4data-v-Ci
   #:sgp4data-v-Di
   #:sgp4data-v-deldt-Earth-zonal-harmonics
   #:sgp4data-deep-space?
   #:sgp4data-resonance
   #:sgp4data-drop-terms?
   #:sgp4data-v-Xi-lunar
   #:sgp4data-v-Zi-lunar
   #:sgp4data-v-Zij-lunar
   #:sgp4data-v-deldt-lunar
   #:sgp4data-v-Xi-solar
   #:sgp4data-v-Zi-solar
   #:sgp4data-v-Zij-solar
   #:sgp4data-v-deldt-solar
   #:sgp4data-v-deltai
   #:sgp4data-v-Dlmpq

   ;; EARTH GRAVITY MODELS
   #:make-v-grav-wgs-72-low-precision
   #:make-v-grav-wgs-72
   #:make-v-grav-wgs-84

   ;; SGP4 INITIALIZATION
   #:init-sgp4

   ;; SGP4 PROPAGATION
   #:propagate-sgp4

   ;; STATE VECTOR UNIT CONVERSIONS
   #:*km-per-unit*
   #:*km/s-per-unit*
   #:*m-per-unit*
   #:*m/s-per-unit*
   #:with-v-grav-units
   #:->v-state-km-s
   #:->v-location-km
   #:->v-velocity-km/s
   #:->v-state-m-s
   #:->v-location-m
   #:->v-velocity-m/s

   ;; CONDITIONS
   #:altitude-lt-0-error
   #:Keplers-equation-convergence-error
   ))
