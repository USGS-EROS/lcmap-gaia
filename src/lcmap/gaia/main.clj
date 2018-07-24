(ns lcmap.gaia.main
  (:gen-class)
  (:require [clojure.tools.logging :as log]
            [lcmap.gaia.config     :as config]
            [lcmap.gaia.util       :as util]))

;; ccdc results are the number of change segments
;; detected over a 100x100 pixel area (chip)
;; There can be multiple changes over time 
;; detected.
;; each pixel will have at least 1 entry
;; if no changes are detected, a result is 
;; returned with an sday and eday == 0


(defmulti gen-product
  (fn [x y product queryday] (keyword product)))

(defmethod gen-product :default
  [x y product queryday]
  (log/infof "%s is not a valid product" product)
  nil)

(defmethod gen-product :primary-landcover
  [x y product queryday])

(defmethod gen-product :secondary-landcover
  [x y product queryday]
  nil)

(defmethod gen-product :primary-confidence
  [x y product queryday]
  nil)

(defmethod gen-product :secondary-confidence
  [x y product queryday]
  nil)

(defmethod gen-product :time-of-change
  [x y product queryday]
  ;; Output the numeric day of year, in the year which a break occurs in.
  ;; Pseudo Code
  ;; queryday = date(“july”, 1, 1990).toordinal()
  ;; if queryday.year == ccdmodel.breakdate.year and ccdmodel.change_prob == 1:
  ;; 	return dayofyear(ccdmodel.breakdate)
  ;; else:
  ;; 	return 0

  ;; queryday -> ordinal format
  (let [query-datetime (util/to-javatime queryday)
        query-year     (util/javatime-year query-datetime)
        break-datetime (util/ordinal-to-javatime 724394)
        break-year     (util/javatime-year break-datetime)

])



)

(defn -main
  ([]
   (println "Welcome to lcmap-gaia\nHere are your options:")
   (System/exit 0))
  ([x y product queryday]
  (log/infof "output for %s %s %s" x y product)
  (System/exit 0)))

(defn foo
  [i]
  (* i 2))
