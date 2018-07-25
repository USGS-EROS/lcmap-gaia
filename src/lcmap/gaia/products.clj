(ns lcmap.gaia.products
  (:require [clojure.tools.logging :as log]
            [lcmap.gaia.util       :as util]))

(defn time-of-change
  [model queryday]
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
        break-day      (get model "bday")
        break-datetime (util/ordinal-to-javatime break-day)
        break-year     (util/javatime-year break-datetime)
        year-match     (= query-year break-year)
        change-prob    (= 1 (get model "chprob"))
        day-of-year    (util/javatime-day-of-year break-day)]
    (if (= true year-match change-prob)
      day-of-year
      0)))
