(ns lcmap.gaia.products-test
  (:require [clojure.test :refer :all]
            [lcmap.gaia.products :as products]
            [lcmap.gaia.file     :as file]
            [lcmap.gaia.util     :as util]))

(def chipdata (file/read-json "resources/y3161805_x-2115585_nodates.json"))
(def pixel_segments (util/pixel-groups chipdata))
(def querydate "2006-07-01")
(def first_pixel (first pixel_segments))
(def pixel_models (last first_pixel))
(def response_set (set [:pixelx :pixely :val]))

(deftest time-of-change-single-model-test
  (let [result (products/time-of-change (first pixel_models) querydate 100 -100)]
    (is (= (set (keys result))  response_set))))

(deftest time-of-change-chip-level-test
  (let [results (map #(products/time-of-change (first %) (last %) querydate) pixel_segments)
        first_result (first results)]
    (is (= (count results) 10000))
    (is (= (set (keys first_result)) response_set))))

(deftest time-since-change-single-model-test
  (let [result (products/time-since-change (first pixel_models) querydate 100 -100)]
    (is (= (set (keys result)) response_set))))

(deftest time-since-change-chip-level-test
  (let [results (map #(products/time-since-change (first %) (last %) querydate) pixel_segments)
        non_nils (filter (fn [i] (some? (:val i))) results)]
    (is (= (count non_nils) 763))
    (is (= (count results) 10000))))

(deftest magnitude-of-change-single-model-test
  (let [result (products/magnitude-of-change (first pixel_models) querydate 100 -100)]
    (is (= (set (keys result)) response_set))))

(deftest magnitude-of-change-chip-level-test
  (let [results (map #(products/magnitude-of-change (first %) (last %) querydate) pixel_segments)
        non_nils (filter (fn [i] (some? (:val i))) results)]
    (is (= (count non_nils) 763))
    (is (= (count results) 10000))))
