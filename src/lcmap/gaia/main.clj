(ns lcmap.gaia.main
  (:gen-class)
  (:require [clojure.tools.logging :as log]
            [lcmap.gaia.config     :as config]
            [lcmap.gaia.file       :as file]
            [lcmap.gaia.util       :as util]))

;; ccdc results are the number of change segments
;; detected over a 100x100 pixel area (chip)
;; There can be multiple changes over time 
;; detected.
;; each pixel will have at least 1 entry
;; if no changes are detected, a result is 
;; returned with an sday and eday == 0


(defmulti gen-product
  (fn [infile product queryday] (keyword product)))

(defmethod gen-product :default
  [infile product queryday]
  (log/infof "%s is not a valid product" product)
  nil)

(defmethod gen-product :primary-landcover
  [infile product queryday])

(defmethod gen-product :secondary-landcover
  [infile product queryday]
  nil)

(defmethod gen-product :primary-confidence
  [infile product queryday]
  nil)

(defmethod gen-product :secondary-confidence
  [infile product queryday]
  nil)

(defmethod gen-product :time-of-change
  [infile product queryday]
  (let [input (file/read-json infile)
        output (map products/time-of-change input)]
      (first output)))

(defn -main
  ([]
   (println "Welcome to lcmap-gaia\nHere are your options:")
   (System/exit 0))
  ([infile product queryday]
  (log/infof "output for %s %s %s" infile product)
  (System/exit 0)))

(defn foo
  [i]
  (* i 2))
