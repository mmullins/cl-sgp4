;;;; -*- Mode: Lisp -*-

(in-package :cl-user)

(defpackage cl-sgp4-system
  (:use :common-lisp :asdf))

(in-package :cl-sgp4-system)

(defsystem cl-sgp4
    :description "SGP4 satellite orbit propagation model"
    :version "0.1"
    :author "Mayes Mullins <mmullins@mullinsenterprises.ca>"
    :licence "MIT"
    :depends-on ()
    :components ((:file "packages")
                 (:file "astrometry" :depends-on ("packages"))
                 (:file "tle" :depends-on ("packages" "astrometry"))
                 (:file "sgp4" :depends-on ("packages" "astrometry"))))
