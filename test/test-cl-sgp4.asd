;;;; -*- Mode: Lisp -*-

(in-package :cl-user)

(defpackage test-cl-sgp4-system
  (:use :common-lisp :asdf))

(in-package :test-cl-sgp4-system)

(defsystem test-cl-sgp4
    :description "Validation functions and sample code for cl-sgp4 package"
    :version "0.1"
    :author "Mayes Mullins <mmullins@mullinsenterprises.ca>"
    :licence "MIT"
    :depends-on (:cl-sgp4)
    :components ((:file "packages")
                 (:file "test-cl-sgp4" :depends-on ("packages"))))
