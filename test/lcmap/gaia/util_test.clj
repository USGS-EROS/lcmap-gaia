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

(deftest to-ordinal-test
  (is (= 726659 (util/to-ordinal date_string))))

(deftest coll-groups-test
  (let [arbitrary_coll [{:a 1 :b 2 :c 3} {:a 1 :b 3 :c 4} {:a 2 :b 5 :c 5}]
        response (util/coll-groups arbitrary_coll [:a])]
    (is (= 2 (count response)))
    (is (= 2 (count (get response {:a 1}))))))

(deftest flatten-vals-test
  (let [coll '({3116865 [{:pixely 3116865, :pixelx -2115585, :val 8} {:pixely 3116865, :pixelx -2115555, :val 9}]})
        resp (util/flatten-vals coll :val)]
    (is (= '(8 9) resp))))

(deftest matching-keys-return-collection-test
  (let [map_a {:foo true :bar false}
        map_b {:foo false :bar true}
        expected_return [map_a map_b]]
    (is (= (util/matching-keys map_a map_b :foo :bar true) expected_return))))

(deftest matching-keys-return-first-map-test
  (let [map_a {:foo true :bar false}
        map_b {:foo false :bar true}
        map_coll [map_a]]
    (is (= (util/matching-keys map_coll map_b :foo :bar true) map_coll))))

(deftest matching-keys-return-last-map-test
  (let [map_a {:foo true :bar false}
        map_b {:foo false :bar true}]
    (is (= (util/matching-keys map_a map_b :foo :bar false) map_b))))

(deftest mean_test
  (let [coll [4 2 6 88 7]]
    (is (= (float 21.4) (util/mean coll)))))





