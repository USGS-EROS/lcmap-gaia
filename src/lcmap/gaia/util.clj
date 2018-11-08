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
  (when (not (nil? datestring))
    (-> datestring (to-javatime) (javatime-to-ordinal))))

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

(defn variable-juxt
  "Return a juxt function that returns the values for the specified map keys"
  [mapkeys]
  (apply juxt (map (fn [i] #(get % i)) mapkeys)))

(defn merge-maps-by-keys
  "Merge two lists of hash-maps, joining the list members by the key values
   specified in the mapkeys argument"
  [maplist1 maplist2 mapkeys]
  (let [conc_lists (concat maplist1 maplist2)
        grouped_lists (group-by (variable-juxt mapkeys) conc_lists)]
    (map #(merge (first %) (last %)) (vals grouped_lists))))

(defn ismap?
  "Returns boolean true / false if input is a hash-map "
  [input]
  (= (type input) clojure.lang.PersistentArrayMap))

(defn matching-keys
  "Return a collection of the map arguments if the key values equal the 
   desired match value, else return map_b. Used in a call to reduce for
   identifying desired maps in a collection"
  [map_a map_b key_a key_b match_value]
  (if (ismap? map_a)
      (if (= match_value (key_a map_a) (key_b map_b))
          (do [map_a map_b])
          (do map_b))       
      (do map_a)))

(defn sort-by-key [coll key] (sort-by (fn [i] (get i key)) coll))



