(ns lcmap.gaia.util
  (:require [java-time :as jt]))

(def gregorian_day_one (jt/local-date 0001 1))

(defn to-javatime
  "Convert a date string to java-time"
  [datestring]
  (jt/local-date "MM-yyyy-dd" datestring))

(defn javatime-year
  "Return year on a java-time object"
  [javatime]
  (if (nil? javatime)
    nil
    (jt/as javatime :year)))

(defn javatime-day-of-year
  "Return day-of-year on a java-time object"
  [javatime]
  (if (nil? javatime)
    nil
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
