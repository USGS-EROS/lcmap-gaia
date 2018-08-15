(ns lcmap.gaia.products
  (:require [clojure.tools.logging :as log]
            [clojure.string        :as string]
            [clojure.math.numeric-tower :as math]
            [lcmap.gaia.util       :as util]))

(defn product-name
  [model product fmt]
  (let [chipx (get model "chipx")
        chipy (get model "chipy")
        name  (string/join "_" [product chipx chipy])]
    (str name "." fmt)))

(defn time-of-change
  "Return numeric day of year in which a break occurs"
  ([model query-day x y]
   (let [change-prob (get model "chprob")
         break-day   (get model "bday")
         query-year  (-> query-day (util/to-javatime) (util/javatime-year))
         break-year  (-> break-day (util/ordinal-to-javatime) (util/javatime-year))
         response    #(hash-map :pixelx x :pixely y :val %)]
     (if (= true (= query-year break-year) (= 1.0 change-prob))
       (-> break-day (util/ordinal-to-javatime) (util/javatime-day-of-year) response) 
       (-> 0 response))))
  ([pixel_map pixel_models query-day]
   (let [values (map #(time-of-change % query-day (get pixel_map "pixelx") (get pixel_map "pixely")) pixel_models)]
     (last (sort-by :val values)))))

(defn time-since-change
  "Return cumulative distance to previous break"
  ([model query-day x y]
   (let [change-prob (get model "chprob")
         break-day   (get model "bday")
         query-ord   (-> query-day (util/to-javatime) (util/javatime-to-ordinal))
         distance    (if (= 1.0 change-prob) (- query-ord break-day) nil)]
     (hash-map :pixelx x :pixely y :val distance)))
  ([pixel_map pixel_models query-day]
   (let [values (map #(time-since-change % query-day (get pixel_map "pixelx") (get pixel_map "pixely")) pixel_models)]
     (first (filter (fn [i] (some? (:val i))) (sort-by :val values))))))

(defn magnitude-of-change
  "Return severity of spectral shift"
  ([model query-day x y]
   (let [change-prob (get model "chprob")
         break-day   (get model "bday")
         query-year  (-> query-day (util/to-javatime) (util/javatime-year))
         break-year  (-> break-day (util/ordinal-to-javatime) (util/javatime-year))
         magnitudes  [(get model "grmag") (get model "remag") (get model "nimag") (get model "s1mag") (get model "s2mag")]
         euc-norm    (math/sqrt (reduce + (map #(math/expt % 2) magnitudes)))
         response    #(hash-map :pixelx x :pixely y :val %)]
     (if (= true (= query-year break-year) (= 1.0 change-prob))
       (-> euc-norm (response))
       (-> 0 (response)))))
  ([pixel_map pixel_models query-day]
   (let [values (map #(magnitude-of-change % query-day (get pixel_map "pixelx") (get pixel_map "pixely")) pixel_models)]
     (last (sort-by :val values)))))

(defn length-of-segment
  "Return length of change segment in days"
  ([model query-day x y]
   (let [query-ord (-> query-day (util/to-javatime) (util/javatime-to-ordinal))
         startends [(- query-ord (get model "sday")) (- query-ord (get model "eday"))]
         positives (filter (fn [i] (> i 0)) startends)
         minimum   (if (= 0 (count positives)) 0 (apply min positives))]
     (hash-map :pixelx x :pixely y :val minimum)))
  ([pixel_map pixel_models query-day]
   (let [values (map #(length-of-segment % query-day (get pixel_map "pixelx") (get pixel_map "pixely")) pixel_models)]
     (first (sort-by :val values)))))
