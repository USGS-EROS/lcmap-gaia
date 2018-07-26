(ns lcmap.gaia.util-test
  (:require [clojure.test :refer :all]
            [lcmap.gaia.util :as util]))

(def date_string "1990-07-12")
(def date_time_string "2018-06-25T20:52:43.508205")

(deftest to-javatime-nil-test
  (is (= (util/to-javatime "not a time") nil)))

(deftest to-javatime-date-test
  (is (= (class (util/to-javatime date_string))
         java.time.LocalDate)))

(deftest to-javatime-date-time-test
  (is (= (class (util/to-javatime date_time_string))
         java.time.LocalDateTime)))

(deftest javatime-year-test
  (let [jt (util/to-javatime date_string)]
    (is (= (util/javatime-year jt) 1990))))

(deftest javatime-day-of-year-test
  (let [jt (util/to-javatime date_string)]
    (is (= (util/javatime-day-of-year jt) 193))))

(deftest ordinal-to-javatime-test
  (is (= (class (util/ordinal-to-javatime 726583)) 
         java.time.LocalDate)))

(deftest javatime-to-ordinal-test
  (let [jt (util/to-javatime date_string)]
    (is (= (util/javatime-to-ordinal jt)
           726659))))

