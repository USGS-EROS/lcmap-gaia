(ns lcmap.gaia.change-products
  (:gen-class)
  (:require [clojure.core.async    :as async]
            [clojure.math.numeric-tower :as math]
            [clojure.math.combinatorics :as combo]
            [clojure.string        :as string]
            [clojure.tools.logging :as log]
            [java-time             :as jt]
            [lcmap.gaia.config     :refer [config]]
            [lcmap.gaia.file       :as file]
            [lcmap.gaia.gdal       :as gdal]
            [lcmap.gaia.product-specs :as product-specs]
            [lcmap.gaia.storage    :as storage]
            [lcmap.gaia.util       :as util]
            [cheshire.core         :as json]))

(defn product-exception-handler
  [exception product_name]
  (let [msg (format "problem calculating %s: %s" product_name (.getMessage exception))]
    (log/error msg)
    (throw (ex-info msg {:type "data-generation-error" :message msg} (.getCause exception)))))

(defn time-of-change
  "Return numeric day of year in which a break occurs"
  ([model date]
   (try
     (let [change-prob (get model "chprob")
           break-day   (util/to-javatime (get model "bday")) 
           query-year  (-> date (util/ordinal-to-javatime) (util/javatime-year))
           break-year  (-> break-day (util/javatime-year))]
       (if (and (= query-year break-year) (= 1.0 change-prob))
         (util/javatime-day-of-year break-day)
         0))
     (catch Exception e
       (product-exception-handler e "time-of-change"))))
  ([pxpy segments date]
   (let [valid (product-specs/segments-valid? segments)
         values (map #(time-of-change % date) segments)]
     (if valid
       (last (sort values))
       0))))

(defn time-since-change
  "Return cumulative distance to previous break"
  ([model date]
   (try
     (let [change-prob (= 1.0 (get model "chprob")) 
           break-day   (util/to-ordinal (get model "bday"))
           day-diff    (- date break-day)]
       (if (and change-prob (>= day-diff 0)) 
                         day-diff 
                         nil))
     (catch Exception e
       (product-exception-handler e "time-since-change"))))
  ([pxpy segments date]
   (let [valid (product-specs/segments-valid? segments)
         values (map #(time-since-change % date) segments)
         not_nil (filter number? values)]
     (if (or (empty? not_nil) (not valid))
       0
       (first (sort not_nil))))))

(defn magnitude-of-change
  "Return severity of spectral shift"
  ([model date]
   (try
     (let [change-prob (= 1.0 (get model "chprob")) 
           query-year  (-> date (util/ordinal-to-javatime) (util/javatime-year))
           break-year  (-> (get model "bday") (util/to-javatime) (util/javatime-year))
           magnitudes  [(get model "grmag") (get model "remag") (get model "nimag") (get model "s1mag") (get model "s2mag")]
           euc-norm    (math/sqrt (reduce + (map #(math/expt % 2) magnitudes)))]
       (if (and (= query-year break-year) change-prob)
         euc-norm
         0))
     (catch Exception e
       (product-exception-handler e "magnitude-of-change"))))
  ([pxpy segments date]
   (let [valid (product-specs/segments-valid? segments)
         values (map #(magnitude-of-change % date) segments)]
     (if valid
       (last (sort values))
       0))))

(defn length-of-segment
  "Return length of change segment in days"
  ([model date]
   (try
     (let [fill (- date (util/to-ordinal (:stability_begin config)))
           start-day (util/to-ordinal (get model "sday"))
           end-day   (util/to-ordinal (get model "eday"))
           diff      (if (> date end-day) (- date end-day) (- date start-day))]
       (if (and (<= 0 diff) (< diff fill)) 
         diff 
         fill))
     (catch Exception e
       (product-exception-handler e "length-of-segment"))))
  ([pxpy segments date]
   (let [fill   (- date (util/to-ordinal (:stability_begin config)))
         valid (product-specs/segments-valid? segments)
         values (map #(length-of-segment % date) segments)]
     (if valid
       (first (sort values))
       fill))))

(defn curve-fit
  "Return Curve QA for point in time"
  ([model date]
   (try
     (let [curve-qa  (get model "curqa")
           start-day (util/to-ordinal (get model "sday"))
           end-day   (util/to-ordinal (get model "eday")) ]
       (if (<= start-day date end-day)
         curve-qa
         0))
     (catch Exception e
       (product-exception-handler e "curve-fit"))))
  ([pxpy segments date]
   (let [valid (product-specs/segments-valid? segments)
         values (map #(curve-fit % date) segments)]
     (if valid
       (last (sort values))
       0))))

(defn products
  [pxpy segments date]
  (let [[px py] pxpy]
    (hash-map :px px :py py :date date
              :values {:curve-fit (curve-fit pxpy segments date)
                       :length-of-segment (length-of-segment pxpy segments date)
                       :magnitude-of-change (magnitude-of-change pxpy segments date)
                       :time-since-change (time-since-change pxpy segments date)
                       :time-of-change (time-of-change pxpy segments date)})))

(defn start-consumers
  [in-chan out-chan operation]
  (let [chunk_size  (:product_instance_count config)]
    (dotimes [_ chunk_size]
      (async/thread
        (while true
          (let [input (async/<!! in-chan)
                result (operation input)]
            (async/>!! out-chan result)))))))

(defn generate
  [{dates :dates cx :cx cy :cy tile :tile :as all}]
  (try
    (let [grouped_segments (storage/grouped-segments cx cy)
          in-chan          (async/chan)
          out-chan         (async/chan)
          operation        #(products (last %) (get grouped_segments (last %)) (first %))
          consumers        (start-consumers in-chan out-chan operation)
          ordinal_dates    (doall (map util/to-ordinal dates))
          pixels           (keys grouped_segments)
          output_fn        (fn [i] (let [result (async/<!! out-chan)] result))
          dates_pixels     (combo/cartesian-product ordinal_dates pixels)]

      (async/go
        (doseq [dp dates_pixels]
          (async/>!! in-chan dp)))

      (let [results (doall (map output_fn dates_pixels))
            date_groups (group-by :date results)] ; group results by (:date )
        (doseq [dg date_groups
                :let [date (util/to-yyyy-mm-dd (first dg))
                      path (storage/ppath "change" cx cy tile date)
                      values (util/flatten-values (second dg))]]
          (util/with-retry (storage/put_json path values))))

      {:products "change" :cx cx :cy cy :dates dates})
    (catch Exception e
      (let [msg (format "problem generating change products for %s: %s" all (.getMessage e))]
        (log/error msg)
        (throw (ex-info msg {:type "data-generation-error" :message msg} (.getCause e)))))))
