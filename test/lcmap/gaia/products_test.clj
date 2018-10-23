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

(def segments (file/read-json "resources/cx-2115585_cy3119805_segment.json"))
(def predictions (file/read-json "resources/cx-2115585_cy3119805_prediction.json"))

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
    (is (= (count gt_zero) 387))
    (is (= (count results) 10000))))

(deftest length-of-segment-single-model-test
  (let [result (products/length-of-segment (first pixel_models) querydate 100 -100)]
    (is (= (set (keys result)) response_set))))

(deftest length-of-segment-chip-level-test
  (let [results (map #(products/length-of-segment (first %) (last %) querydate) tr/pixel_segments)
        gt_zero (filter (fn [i] (> (:val i) 0)) results)]
    (is (= (count gt_zero) 5633))
    (is (= (count results) 10000))))

(deftest curve-fit-single-model-test
  (let [result (products/curve-fit (first pixel_models) querydate 100 -100)]
    (is (= (set (keys result)) response_set))))

(deftest curve-fit-chip-level-test
  (let [results (map #(products/curve-fit (first %) (last %) querydate) tr/pixel_segments)
        gt_zero (filter (fn [i] (> (:val i) 0)) results)]
    (is (= (count gt_zero) 9989))
    (is (= (count results) 10000))))

;; (deftest data-test
;;   ;; FIXME !!
;;   (let [values (products/data chip_data chip_data "time-since-change" querydate)] 
;;     (is (= (count values) 10000))))

(deftest ismap?-affirmative-test
  (is (products/ismap? {:a "map"})))

(deftest ismap?-negative-test
  (is (not (products/ismap? :not_a_map))))

(deftest matching-keys-return-collection-test
  (let [map_a {:foo true :bar false}
        map_b {:foo false :bar true}
        expected_return [map_a map_b]]
    (is (= (products/matching-keys map_a map_b :foo :bar true) expected_return))))

(deftest matching-keys-return-first-map-test
  (let [map_a {:foo true :bar false}
        map_b {:foo false :bar true}
        map_coll [map_a]]
    (is (= (products/matching-keys map_coll map_b :foo :bar true) map_coll))))

(deftest matching-keys-return-last-map-test
  (let [map_a {:foo true :bar false}
        map_b {:foo false :bar true}]
    (is (= (products/matching-keys map_a map_b :foo :bar false) map_b))))

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
