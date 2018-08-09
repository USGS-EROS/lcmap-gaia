(ns lcmap.gaia.products
  (:require [clojure.tools.logging :as log]
            [clojure.string        :as string]
            [lcmap.gaia.util       :as util]))

(defn product-name
  [model product fmt]
  (let [chipx (get model "chipx")
        chipy (get model "chipy")
        name  (string/join "_" [product chipx chipy])]
    (str name "." fmt)))

(defn time-of-change
  ([model query-day x y]
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
        break-year  (-> break-day (util/ordinal-to-javatime) (util/javatime-year))
        response    #(hash-map :pixelx x :pixely y :toc %)]
    (if (= true (= query-year break-year) (= 1.0 change-prob))
      (-> break-day (util/ordinal-to-javatime) (util/javatime-day-of-year) response) 
      (-> 0 response))))
  ([pixel_map pixel_models query-day]
   (let [values (map #(time-of-change % query-day (:pixelx pixel_map) (:pixely pixel_map)) pixel_models)]
     (last (sort-by :toc values)))))


