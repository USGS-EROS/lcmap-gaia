(ns lcmap.gaia.products-test
  (:require [clojure.test :refer :all]
            [lcmap.gaia.products :as products]
            [lcmap.gaia.file     :as file]
            [lcmap.gaia.util     :as util]
            [lcmap.gaia.test-resources :as tr]))

(def querydate "2006-07-01")
(def chip_data tr/chip_data)
(def first_pixel (first tr/pixel_segments))
(def pixel_models (last first_pixel))
(def response_set (set [:pixelx :pixely :val]))

(deftest time-of-change-single-model-test
  (let [result (products/time-of-change (first pixel_models) querydate 100 -100)]
    (is (= (set (keys result))  response_set))))

(deftest time-of-change-chip-level-test
  (let [results (map #(products/time-of-change (first %) (last %) querydate) tr/pixel_segments)
        first_result (first results)]
    (is (= (count results) 10000))
    (is (= (set (keys first_result)) response_set))))

(deftest time-since-change-single-model-test
  (let [result (products/time-since-change (first pixel_models) querydate 100 -100)]
    (is (= (set (keys result)) response_set))))

(deftest time-since-change-chip-level-test
  (let [results (map #(products/time-since-change (first %) (last %) querydate) tr/pixel_segments)
        non_nils (filter (fn [i] (some? (:val i))) results)]
    ;(is (= (count non_nils) 763))
    (is (= (count results) 10000))))

(deftest magnitude-of-change-single-model-test
  (let [result (products/magnitude-of-change (first pixel_models) querydate 100 -100)]
    (is (= (set (keys result)) response_set))))

(deftest magnitude-of-change-chip-level-test
  (let [results (map #(products/magnitude-of-change (first %) (last %) querydate) tr/pixel_segments)
        gt_zero (filter (fn [i] (> (:val i) 0)) results)]
    (is (= (count gt_zero) 54))
    (is (= (count results) 10000))))

(deftest length-of-segment-single-model-test
  (let [result (products/length-of-segment (first pixel_models) querydate 100 -100)]
    (is (= (set (keys result)) response_set))))

(deftest length-of-segment-chip-level-test
  (let [results (map #(products/length-of-segment (first %) (last %) querydate) tr/pixel_segments)
        gt_zero (filter (fn [i] (> (:val i) 0)) results)]
    (is (= (count gt_zero) 9942))
    (is (= (count results) 10000))))

(deftest curve-fit-single-model-test
  (let [result (products/curve-fit (first pixel_models) querydate 100 -100)]
    (is (= (set (keys result)) response_set))))

(deftest curve-fit-chip-level-test
  (let [results (map #(products/curve-fit (first %) (last %) querydate) tr/pixel_segments)
        gt_zero (filter (fn [i] (> (:val i) 0)) results)]
    (is (= (count gt_zero) 9977))
    (is (= (count results) 10000))))

(deftest data-test
  (let [values (products/data chip_data "time-since-change" querydate)] 
    (is (= (count values) 10000))))
