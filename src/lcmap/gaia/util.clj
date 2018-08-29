(ns lcmap.gaia.util
  (:require [java-time :as jt]
            [cheshire.core :as json]
            [org.httpkit.client :as http]
            [lcmap.gaia.config :refer [config]]))

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

(defn coll-groups
  "Group collection of hash maps by shared keys values"
  [coll keys]
  (group-by #(select-keys % keys) coll))

(defn pixel-groups
  [coll]
  (coll-groups coll ["pixelx" "pixely"]))

(defn flatten-vals
  "Flatten the values for a collection of hash-maps"
  [coll mapkey]
  (let [coll_vals (map (fn [i] (vals i)) coll)
        vals_flat (flatten coll_vals)]
    (map mapkey vals_flat)))

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
  []
  (let [grid_resource (str (:chipmunk_host config) "/grid")
        grid_response (http/get grid_resource)
        json_body (first (json/parse-string (:body @grid_response)))]
    (get json_body "proj")))
