(ns lcmap.gaia.util
  (:gen-class)
  (:require [java-time :as jt]
            [cheshire.core :as json]
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

(defn to-ordinal
  "Convert ISO8601 date string to an ordinal value"
  [datestring]
  (-> datestring (to-javatime) (javatime-to-ordinal)))

(defn coll-groups
  "Group collection of hash maps by shared keys values"
  [coll keys]
  (group-by #(select-keys % keys) coll))

(defn pixel-groups
  [coll]
  (coll-groups coll ["px" "py"]))

(defn flatten-vals
  "Flatten the values for a collection of hash-maps"
  [coll mapkey]
  (let [coll_vals (map (fn [i] (vals i)) coll)
        vals_flat (flatten coll_vals)]
    (map mapkey vals_flat)))

