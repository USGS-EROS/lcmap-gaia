(ns lcmap.gaia.util
  (:gen-class)
  (:require [again.core            :as again]
            [cheshire.core         :as json]
            [clojure.string        :as string]
            [clojure.tools.logging :as log]
            [environ.core          :as environ]
            [java-time             :as jt]
            [org.httpkit.client    :as http]
            [lcmap.gaia.file       :as file]
            [lcmap.gaia.config     :refer [config]]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(def gregorian_day_one (jt/local-date 0001 1))
(def date_pattern (re-pattern #"[0-9]{4}-[0-9]{2}-[0-9]{2}"))
(def date_time_pattern (re-pattern #"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}\.[0-9]{6}"))

(defmulti to-javatime
  (fn [datestring]
    (cond (re-matches date_pattern datestring) :year-month-day
          (re-matches date_time_pattern datestring) :iso8601)))

(defmethod to-javatime :default
  [datestring] nil)

(defmethod to-javatime :year-month-day
  [datestring]
  (jt/local-date datestring))

(defmethod to-javatime :iso8601
  [datestring]
  (jt/local-date-time datestring))

(defn javatime-year
  "Return year on a java-time object"
  [javatime]
  (when (not (nil? javatime))
    (jt/as javatime :year)))

(defn javatime-day-of-year
  "Return day-of-year on a java-time object"
  [javatime]
  (when (not (nil? javatime))
    (jt/as javatime :day-of-year)))

(defn ordinal-to-javatime
  "Convert an ordinal day on the Gregorian Calendar
  to java-time"
  [^long ordinal]
  (let [days_from_zero (- ordinal 1)]
    (jt/plus gregorian_day_one (jt/days days_from_zero))))

(defn javatime-to-ordinal
  "Convert java-time to ordinal value"
  [javatime]
  (jt/time-between gregorian_day_one javatime :days))

(defn to-ordinal
  "Convert ISO8601 date string to an ordinal value"
  [datestring]
  (when (not (nil? datestring))
    (-> datestring (to-javatime) (javatime-to-ordinal))))

(defn to-yyyy-mm-dd
  "Convert ordinal value to a YYYY-MM-DD formatted string"
  [^long ordinaldate]
  (-> ordinaldate inc ordinal-to-javatime str))

(defn todays-date
  "Return todays date as a string"
  []
  (str (jt/local-date)))

(defn coll-groups
  "Group collection of hash maps by shared keys values"
  [coll keys]
  (group-by #(select-keys % keys) coll))

(defn flatten-vals
  "Flatten the values for a collection of hash-maps"
  [coll mapkey]
  (let [vals_flat (-> coll (#(map vals %)) (flatten))]
    (map mapkey vals_flat)))

(defn variable-juxt
  "Return a juxt function that returns the values for the specified map keys"
  [mapkeys]
  (apply juxt (map (fn [i] #(get % i)) mapkeys)))

(defn matching-keys
  "Return a collection of the map arguments if the key values equal the 
   desired match value, else return map_b. Used in a call to reduce for
   identifying desired maps in a collection"
  [map_a map_b key_a key_b match_value]
  (if (map? map_a)
      (if (= match_value (key_a map_a) (key_b map_b))
          (do [map_a map_b])
          (do map_b))       
      (do map_a)))

(defn sort-by-key [coll key] (sort-by (fn [i] (get i key)) coll))

(defn subtract_year
  "Return an ordinal date for one year prior to provided value"
  [ordinal_date]
  (let [jt_date (ordinal-to-javatime ordinal_date)
        minus_year (jt/minus jt_date (jt/years 1))   ; if jt_date is July 1, minus_year will be June 30th
        plus_day   (jt/plus minus_year (jt/days 1))] ; add day to get to July 1
    (javatime-to-ordinal plus_day)))

(defn concat_ints
  [& args]
  (let [str_ints (map str args)
        concated (string/join str_ints)]
    (-> concated (bigdec) (int))))

(defn scale-value
  "Return scaling of probability into integer, with a min value of 1"
  ([^double value ^long factor]
    (let [_prob (* value factor)]
      (if (< _prob 1)
        1
        (int _prob))))
  ([value]
   (scale-value value 100)))

(defn mean 
  "Returns the mathematical mean value for a collection of numbers"
  [coll]
  (let [sum (float (apply + coll)) 
        cnt (int (count coll))]
    (if (pos? cnt)
      (float (/ sum cnt)) 
      0)))

;; add-usr-path and amend-usr-path blatantly ripped off from the
;; USGS-EROS/lcmap-chipmunk project on GitHub, created by
;; Jon Morton https://github.com/jmorton
;; 
(defn add-usr-path
  ""
  [& paths]
  (let [field (.getDeclaredField ClassLoader "usr_paths")]
    (try (.setAccessible field true)
         (let [original (vec (.get field nil))
               updated  (distinct (concat original paths))]
           (.set field nil (into-array updated)))
         (finally
           (.setAccessible field false)))))

(defn amend-usr-path
  ""
  [more-paths]
  (apply add-usr-path more-paths))

(defn get-projection
  ([]
   (try
     (let [grid_resource (str (:chipmunk-host environ/env) "/grid")
           grid_response (http/get grid_resource)
           json_body (first (json/parse-string (:body @grid_response)))]
       (get json_body "proj"))
     (catch Exception e
       (let [grid_resource (str (:chipmunk-host environ/env) "/grid")
             msg (format "problem in util/get-projection %s: %s" grid_resource (.getMessage e))]
         (log/error msg)
         (throw (ex-info msg {:type "data-request-error" :message msg} (.getCause e)))))))
  ([local]
   (let [grid (first (file/read-json "resources/grid.conus.json"))]
     (get grid "proj"))))

(defn float-string
  "Ensure value is in float format"
  [input]
  (if (number? input)
    (-> input (float) (str))
    (-> input (read-string) (float) (str))))

(defn pixel-groups
  [injson]
  (let [juxt_fn (variable-juxt ["px" "py"])]
    (group-by juxt_fn injson)))

(defn flatten-values
  "Return a flat list of product values given a collection of hash-maps
  for every pixel in a chip, [{:pixely 3159045, :pixelx -2114775, :val 6290},...] ...]"
  [product_value_collection]
  (try
    (let [; group product coll by row
          row_groups (coll-groups product_value_collection [:py]) 
          ; sort row group values by pixelx ascending 
          sort-pixelx-fn (fn [i] (hash-map (:py (first i)) (sort-by :px (last i))))
          sorted-x-vals (map sort-pixelx-fn row_groups)
          ; sort the rows by the pixely key ascending
          sorted-y-rows (sort-by (fn [i] (first (keys i))) > sorted-x-vals)
          ; finally, flatten to a one dimensional list
          flattened (flatten-vals sorted-y-rows :values)]
      flattened)
    (catch Exception e
      (let [msg (format "problem in util/flatten_values - input count: %s first input: %s last input: %s message: %s"
                        (count product_value_collection) (first product_value_collection) (last product_value_collection) (.getMessage e))]
        (log/error msg)
        (throw (ex-info msg {:type "data-generation-error" :message msg} (.getCause e)))))))

(defmacro log-time
  "Evaluates expr and prints the time it took.  Returns the value of
 expr." ; poached from clojure.core
  [expr description]
  `(let [start# (. System (nanoTime))
         ret# ~expr]
     (log/info (str ~description " took: " (/ (double (- (. System (nanoTime)) start#)) 1000000.0) " msecs"))
     ret#))

(defmacro with-retry
  [expr]
  `(again/with-retries (:retry_strategy config) ~expr))

