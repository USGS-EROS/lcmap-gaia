(ns lcmap.gaia.products-test
  (:require [clojure.test :refer :all]
            [lcmap.gaia.products :as products]
            [lcmap.gaia.file     :as file]
            [lcmap.gaia.util     :as util]
            [lcmap.gaia.test-resources :as tr]))


(def first_pixel (first tr/pixel_map))
(def first_segments_predictions (first (vals first_pixel))) 
(def response_set (set [:pixelx :pixely :val]))

(deftest time-of-change-single-model-test
  (let [result (products/time-of-change (first (:segments first_segments_predictions))  tr/querydate 100 -100)]
    (is (= (set (keys result))  response_set))))

(deftest time-of-change-chip-level-test
  (let [results (map #(products/time-of-change (-> % (keys) (first)) (-> % (vals) (first)) tr/querydate) tr/pixel_map)
        first_result (first results)]
    (is (= (count results) 10000))
    (is (= (set (keys first_result)) response_set))))

(deftest time-since-change-single-model-test
  (let [result (products/time-since-change (first (:segments first_segments_predictions)) tr/querydate 100 -100)]
    (is (= (set (keys result)) response_set))))

(deftest time-since-change-chip-level-test
  (let [results (map #(products/time-since-change (-> % (keys) (first)) (-> % (vals) (first)) tr/querydate) tr/pixel_map)
        non_nils (filter (fn [i] (some? (:val i))) results)]
    ;(is (= (count non_nils) 763))
    (is (= (count results) 10000))))

(deftest magnitude-of-change-single-model-test
  (let [result (products/magnitude-of-change (first (:segments first_segments_predictions)) tr/querydate 100 -100)]
    (is (= (set (keys result)) response_set))))

(deftest magnitude-of-change-chip-level-test
  (let [results (map #(products/magnitude-of-change (-> % (keys) (first)) (-> % (vals) (first)) tr/querydate) tr/pixel_map)
        gt_zero (filter (fn [i] (> (:val i) 0)) results)]
    (is (= (count gt_zero) 387))
    (is (= (count results) 10000))))

(deftest length-of-segment-single-model-test
  (let [result (products/length-of-segment (first (:segments first_segments_predictions)) tr/querydate 100 -100)]
    (is (= (set (keys result)) response_set))))

(deftest length-of-segment-chip-level-test
  (let [results (map #(products/length-of-segment (-> % (keys) (first)) (-> % (vals) (first)) tr/querydate) tr/pixel_map)
        gt_zero (filter (fn [i] (> (:val i) 0)) results)]
    (is (= (count gt_zero) 5633))
    (is (= (count results) 10000))))

(deftest curve-fit-single-model-test
  (let [result (products/curve-fit (first (:segments first_segments_predictions)) tr/querydate 100 -100)]
    (is (= (set (keys result)) response_set))))

(deftest curve-fit-chip-level-test
  (let [results (map #(products/curve-fit (-> % (keys) (first)) (-> % (vals) (first)) tr/querydate) tr/pixel_map)
        gt_zero (filter (fn [i] (> (:val i) 0)) results)]
    (is (= (count gt_zero) 9989))
    (is (= (count results) 10000))))

;; (deftest data-test
;;   ;; FIXME !!
;;   (let [values (products/data chip_data chip_data "time-since-change" querydate)] 
;;     (is (= (count values) 10000))))

(deftest falls_between_eday_sday-coll-test
  (let [map_a {:follows_eday true :precedes_sday false}
        map_b {:precedes_sday true :follows_eday false}
        expected [map_a map_b]]
    (is (= (products/falls-between-eday-sday map_a map_b) expected))))

(deftest falls_between_eday_sday-coll-test
  (let [map_a {:follows_eday true :precedes_sday false}
        map_b {:precedes_sday true :follows_eday false}
        expected [map_a map_b]]
    (is (= (products/falls-between-eday-sday map_a map_b) expected))))

(deftest falls_between_eday_sday-map-test
  (let [map_a {:follows_eday false :precedes_sday true}
        map_b {:precedes_sday true :follows_eday false}]
    (is (= (products/falls-between-eday-sday map_a map_b) map_b))))

(deftest falls_between_eday_sday-nonmap-test
  (let [map_a {:follows_eday true :precedes_sday false}
        map_b {:precedes_sday true :follows_eday false}
        map_c {:precedes_sday true :follows_eday false}
        expected [map_a map_b]]
    (is (= (products/falls-between-eday-sday [map_a map_b] map_c) expected))))

(deftest falls_between_bday_sday-coll-test
  (let [map_a {:follows_bday true :precedes_sday false}
        map_b {:precedes_sday true :follows_bday false}
        expected [map_a map_b]]
    (is (= (products/falls-between-bday-sday map_a map_b) expected))))

;; (deftest model_class_test
;;   (let [model (first first_pixel_models)
;;         query_ord (-> "1996-07-01" (util/to-ordinal))
;;         model_rank_0 (products/model-class model query_ord 0)
;;         model_rank_1 (products/model-class model query_ord 1)]
;;     (print model_rank_0)
;;     (is (= 0 (:value model_rank_0)))
;;     (is (false? (:intersects model_rank_0)))
;;     (is (= 5 (:value model_rank_1)))
;;     (is (:precedes_sday model_rank_0))))

