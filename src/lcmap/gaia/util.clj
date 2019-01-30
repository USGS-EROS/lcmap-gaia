(ns lcmap.gaia.util
  (:gen-class)
  (:require [java-time :as jt]
            [cheshire.core :as json]
            [clojure.string :as string]
            [lcmap.gaia.config :refer [config]]
            [environ.core :as environ]
            [org.httpkit.client :as http]
            [lcmap.gaia.file :as file]))

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
  [ordinal]
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
  ([value factor]
    (let [_prob (* value factor)]
      (if (< _prob 1)
        1
        (int _prob))))
  ([value]
   (scale-value value 100)))

(defn mean 
  "Returns the mathematical mean value for a collection of numbers"
  [coll]
  (let [sum (apply + coll)
        count (count coll)]
    (if (pos? count)
      (float (/ sum count)) 
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
   (let [grid_resource (str (:chipmunk-host environ/env) "/grid")
         grid_response (http/get grid_resource)
         json_body (first (json/parse-string (:body @grid_response)))]
     (get json_body "proj")))
  ([local]
   (let [grid (first (file/read-json "resources/grid.conus.json"))]
     (get grid "proj"))))

(defn product-output-name
  ([product x y date suffix]
   (->> [product x y date] (string/join "-") (#(str % suffix))))
  ([product x y date]
   (product-output-name product x y date ".json")))
