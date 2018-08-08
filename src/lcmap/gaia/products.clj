(ns lcmap.gaia.products
  (:require [clojure.tools.logging :as log]
            [lcmap.gaia.util       :as util]))

(defn time-of-change
  [model query-day]
  ;; Output the numeric day of year, in the year which a break occurs in.
  ;; Pseudo Code
  ;; queryday = date(“july”, 1, 1990).toordinal()
  ;; if queryday.year == ccdmodel.breakdate.year and ccdmodel.change_prob == 1:
  ;; 	return dayofyear(ccdmodel.breakdate)
  ;; else:
  ;; 	return 0

  ;; queryday -> ordinal format
  (let [change-prob (get model "chprob")
        break-day   (get model "bday")
        query-year  (-> query-day (util/to-javatime) (util/javatime-year))
        break-year  (-> break-day (util/ordinal-to-javatime) (util/javatime-year))]
    (if (= true (= query-year break-year) (= 1.0 change-prob))
      (-> break-day (util/ordinal-to-javatime) (util/javatime-day-of-year)) 
      0)))
