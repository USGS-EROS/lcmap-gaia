(ns lcmap.gaia.util
  (:require [java-time :as jt]))

(defn to-javatime
  "Convert a date string to java-time"
  [datestring]
  (jt/local-date "MM-yyyy-dd" datestring))

(defn javatime-year
  "Return year on a java-time object"
  [javatime]
  (jt/as javatime :year))

(defn javatime-day-of-year
  "Return day-of-year on a java-time object"
  [javatime]
  (jt/as javatime :day-of-year))

(defn ordinal-to-javatime
  "Convert an ordinal day on the Gregorian Calendar
  to java-time"
  [ordinal]
  (let [dayone_yearone (jt/local-date 0001 1)
        days_from_zero (- ordinal 1)]
    (jt/plus dayone_yearone (jt/days days_from_zero))))

(defn javatime-to-ordinal
  "Convert java-time to ordinal value"
  [javatime]
  ;;  (jt/time-between (jt/local-date 0001 1) (jt/local-date 2017 9 28) :days)
)
