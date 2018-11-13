(ns lcmap.gaia.products-test
  (:require [clojure.test :refer :all]
            [lcmap.gaia.products :as products]
            [lcmap.gaia.file     :as file]
            [lcmap.gaia.util     :as util]
            [lcmap.gaia.config   :refer [config]]
            [lcmap.gaia.test-resources :as tr]))


(def first_pixel (first tr/pixel_map))
(def first_segments_predictions (first (vals first_pixel))) 
(def response_set (set [:pixelx :pixely :val]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;    CHANGE PRODUCT TESTS    ;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(deftest time-of-change-single-model-test
  (let [result (products/time-of-change (first (:segments first_segments_predictions))  tr/query_ord 100 -100)]
    (is (= (set (keys result))  response_set))))

(deftest time-of-change-chip-level-test
  (let [results (map #(products/time-of-change (-> % (keys) (first)) (-> % (vals) (first)) tr/query_ord) tr/pixel_map)
        first_result (first results)]
    (is (= (count results) 10000))
    (is (= (set (keys first_result)) response_set))))

(deftest time-since-change-single-model-test
  (let [result (products/time-since-change (first (:segments first_segments_predictions)) tr/query_ord 100 -100)]
    (is (= (set (keys result)) response_set))))

(deftest time-since-change-chip-level-test
  (let [results (map #(products/time-since-change (-> % (keys) (first)) (-> % (vals) (first)) tr/query_ord) tr/pixel_map)
        non_nils (filter (fn [i] (some? (:val i))) results)]
    ;(is (= (count non_nils) 763))
    (is (= (count results) 10000))))

(deftest magnitude-of-change-single-model-test
  (let [result (products/magnitude-of-change (first (:segments first_segments_predictions)) tr/query_ord 100 -100)]
    (is (= (set (keys result)) response_set))))

(deftest magnitude-of-change-chip-level-test
  (let [results (map #(products/magnitude-of-change (-> % (keys) (first)) (-> % (vals) (first)) tr/query_ord) tr/pixel_map)
        gt_zero (filter (fn [i] (> (:val i) 0)) results)]
    (is (= (count gt_zero) 387))
    (is (= (count results) 10000))))

(deftest length-of-segment-single-model-test
  (let [result (products/length-of-segment (first (:segments first_segments_predictions)) tr/query_ord 100 -100)]
    (is (= (set (keys result)) response_set))))

(deftest length-of-segment-chip-level-test
  (let [results (map #(products/length-of-segment (-> % (keys) (first)) (-> % (vals) (first)) tr/query_ord) tr/pixel_map)
        gt_zero (filter (fn [i] (> (:val i) 0)) results)]
    (is (= (count gt_zero) 5633))
    (is (= (count results) 10000))))

(deftest curve-fit-single-model-test
  (let [result (products/curve-fit (first (:segments first_segments_predictions)) tr/query_ord 100 -100)]
    (is (= (set (keys result)) response_set))))

(deftest curve-fit-chip-level-test
  (let [results (map #(products/curve-fit (-> % (keys) (first)) (-> % (vals) (first)) tr/query_ord) tr/pixel_map)
        gt_zero (filter (fn [i] (> (:val i) 0)) results)]
    (is (= (count gt_zero) 9989))
    (is (= (count results) 10000))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;    CLASSIFICATION PRODUCT TESTS    ;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

(deftest nbr_test
  (let [first_nbr (products/nbr tr/first_segment)
        last_nbr  (products/nbr tr/last_segment)]
    (is (> first_nbr 0.12))
    (is (< first_nbr 0.14))
    (is (> last_nbr  0.023))
    (is (< last_nbr  0.024))))

(deftest get_class_test
  (let [first_class (products/get-class tr/first_probs)
        last_class  (products/get-class tr/last_probs)]
    (is (= first_class 7))
    (is (= last_class 2))))

(deftest first_date_of_class_test
  (let [sorted_predictions (util/sort-by-key tr/first_grouped_predictions :date)]
    (is (= "1995-07-01" (products/first-date-of-class sorted_predictions 7)))
    (is (= "2012-07-01" (products/first-date-of-class sorted_predictions 5)))
    (is (= nil (products/first-date-of-class sorted_predictions 4)))))

(deftest mean_test
  (let [coll [4 2 6 88 7]]
    (is (= (float 21.4) (products/mean coll)))))

(deftest mean_probabilities_test
  (let [preds [{:prob [0 1 2 3 4 5 6 7 8]} {:prob [7 8 9 5 8 7 6 5 5]}]]
    (is (= [3.5 4.5 5.5 4.0 6.0 6.0 6.0 6.0 6.5]
           (products/mean-probabilities preds)))))

(deftest classify_positive_nbr_test
   (let [model (merge tr/first_segment_modded {:probabilities tr/grass_to_forest_probs})
         post_forest_query_date (-> "2001-07-01" (util/to-ordinal))
         pre_forest_query_date (-> "1998-07-01" (util/to-ordinal))
         nbrdiff (float 0.06)]
     ; default value 3 is grass, 4 is tree
     (is (= 4 (products/classify model post_forest_query_date 0 nbrdiff)))
     (is (= 3 (products/classify model pre_forest_query_date 0 nbrdiff)))
     (is (= 3 (products/classify model post_forest_query_date 1 nbrdiff)))
     (is (= 4 (products/classify model pre_forest_query_date 1 nbrdiff)))))

(deftest classify_negative_nbr_test
  (let [model (merge tr/first_segment_modded {:probabilities tr/forest_to_grass_probs})
        post_grass_query_date (-> "2001-07-01" (util/to-ordinal))
        pre_grass_query_date (-> "1998-07-01" (util/to-ordinal))
        nbrdiff (float -0.06)]
    ; default value 3 is grass, 4 is tree
    (is (= 3 (products/classify model post_grass_query_date 0 nbrdiff)))
    (is (= 4 (products/classify model pre_grass_query_date 0 nbrdiff)))
    (is (= 4 (products/classify model post_grass_query_date 1 nbrdiff)))
    (is (= 3 (products/classify model pre_grass_query_date 1 nbrdiff)))))

(deftest classify_else_test
  (let [first_segment (first tr/first_sorted_segments)
        sday (-> first_segment (:sday) (util/to-ordinal))
        nbrdiff (products/nbr first_segment)
        segment_probabilities (filter (fn [i] (= (util/to-ordinal (:sday i)) sday)) tr/first_probabilities)
        sorted_probabilities (util/sort-by-key segment_probabilities :date)
        segment_model (merge first_segment {:probabilities sorted_probabilities})]
    (is (= 7 (products/classify segment_model tr/query_ord 0 nbrdiff)))))

(deftest characterize_segment_test
  (with-redefs [products/nbr (fn [i] 66)
                products/classify (fn [a b c d] 99)]
    (let [segment {:sday "1990-04-27" :eday "2000-06-11" :bday "2000-06-11"}
          query_day (-> "1998-07-01" (util/to-ordinal))
          probabilities [{:sday "1990-04-27" :date "1995-07-01"} 
                         {:sday "2000-07-10" :date "2001-07-01"}]
          characterized (products/characterize_segment segment query_day probabilities 0)]
      (is (= characterized {:intersects true
                            :precedes_sday false
                            :follows_eday false
                            :follows_bday false
                            :btw_eday_bday false
                            :sday (-> "1990-04-27" (util/to-ordinal))
                            :eday (-> "2000-06-11" (util/to-ordinal))
                            :bday (-> "2000-06-11" (util/to-ordinal))
                            :growth true
                            :decline false
                            :probabilities '({:sday "1990-04-27", :date "1995-07-01"})
                            :classification 99})))))

(deftest landcover_test ; first segment -> sday 1982-12-27 bday 2001-10-04 eday 2001-09-10
                        ; last segment -> sday 2001-10-04  bday 2017-09-14 eday 2017-09-14
  (let [segments_probabilities tr/first_segments_predictions]
    ; query date precedes first segment start date and fill_begin is true
    (is (= (:snow (:lc_map config)) (products/landcover segments_probabilities (-> "1980-01-01" (util/to-ordinal)) 0)))
    
    ; query date precedes first segment start date
    (is (= (:lc_insuff (:lc_defaults config))
           (products/landcover segments_probabilities (-> "1980-01-01" (util/to-ordinal)) 0 (merge config {:fill_begin false}))))

    ; query date follows last segment end date and fill_end is true
    (is (= (:water (:lc_map config)) (products/landcover segments_probabilities (-> "2018-01-01" (util/to-ordinal)) 0)))

    ; query date follows last segment end date
    (is (= (:lc_insuff (:lc_defaults config))
           (products/landcover segments_probabilities (-> "2017-10-01" (util/to-ordinal)) 0 (merge config {:fill_end false}))))

    ; query date falls between a segments start and end dates
    (is (= (:water (:lc_map config))
           (products/landcover segments_probabilities (-> "2002-01-01" (util/to-ordinal)) 0)))

    ; query date falls between segments of same landcover classification and fill_samelc config is true
    (is (= (:tree (:lc_map config))
           (products/landcover tr/first_segments_matching_predictions (-> "2001-09-20" (util/to-ordinal)) 0)))

    (let [modded_segments {:segments [tr/first_segment_modded (last (:segments tr/first_segments_predictions))] :predictions (:predictions tr/first_segments_predictions)}]
      ; query date falls between one segments break date and the following segments start date and fill_difflc config is true
      (is (= (:water (:lc_map config)) (products/landcover modded_segments (-> "2001-10-03" (util/to-ordinal)) 0)))
      ; query date falls between a segments end date and breake date and fill_difflc config is true
      (is (= (:snow (:lc_map config)) (products/landcover modded_segments (-> "2001-09-20" (util/to-ordinal)) 0)))
      ; as a last resort return lc_inbtw configuration value
      (is (= (:lc_inbtw config))
          (products/landcover modded_segments (-> "2001-09-20" (util/to-ordinal)) 0 (merge config {:fill_difflc false}))))))

