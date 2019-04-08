(ns lcmap.gaia.products-test
  (:require [clojure.test :refer :all]
            [clojure.string      :as string]
            [lcmap.gaia.products :as products]
            [lcmap.gaia.file     :as file]
            [lcmap.gaia.util     :as util]
            [lcmap.gaia.config   :refer [config]]
            [lcmap.gaia.storage  :as storage]
            [lcmap.gaia.nemo     :as nemo]
            [lcmap.gaia.test-resources :as tr]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;    CHANGE PRODUCT TESTS    ;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def response_set (set [:pixelx :pixely :val]))

(deftest get-prefix-test
  (let [grid "cu"
        date "2001-07-01"
        tile "123456"
        result (products/get-prefix grid date tile)]
    (is (= result "2001/cu/123/456"))))

(deftest map-path-test
  (with-redefs [config {:region "CU" :ccd_ver "C01"}]
    (let [tileid "123456"
          product "time-since-change"
          date "2007-07-01"
          result (products/map-path tileid product date)]
      (is (= (keys result) '(:name :prefix)))
      (is (= (:prefix result) "2007/CU/123/456"))
      (is (string/includes? (:name result) "LCMAP-CU-123456-20070701-"))
      (is (string/includes? (:name result) "-C01-SCLAST.tif")))))

(deftest ppath-test
  (with-redefs [config {:region "CU"}]
    (let [product "TSC"
          x "111111"
          y "222222"
          tile "345678"
          date "2007-07-01"
          result (products/ppath product x y tile date)]
      (is (= result {:name "TSC-111111-222222-2007-07-01.json", :prefix "2007/CU/345/678"})))))

(deftest persist-test
  (with-redefs [config {:region "CU"}
                storage/put_json (fn [a b] true)
                products/data (fn [a b c d] [1 2 3 4])]
    (let [product "TSC"
          cx "111111"
          cy "222222"
          tile "012345"
          query_day "2007-07-01"
          data {:segments [] :predictions []}
          results (products/persist product cx cy tile query_day data)]
      (is (= results {:cx cx :cy cy :date query_day :status "success"})))))

(deftest generation-test
  (with-redefs [nemo/results (fn [cx cy] [{:segments [1 2 3] :predictions [4 5 6]}])
                products/persist (fn [p cx cy tile day data] {:date day :status "fail" :message "bad message"})]
    (let [input {:dates ["2006-07-01"] :cx 111111 :cy 222222 :product "TSC" :tile "012345"}
          result (products/generation input)]
      (is (= result {:failures '({"2006-07-01" "bad message"}), :product "TSC", :cx 111111, :cy 222222, :dates ["2006-07-01"]})))))

(deftest time-of-change-single-model-test
  (let [result (products/time-of-change (first (:segments tr/first_segments_predictions))  tr/query_ord 100 -100)]
    (is (= (set (keys result))  response_set))))

(deftest time-of-change-chip-level-test
  (let [results (map #(products/time-of-change (-> % (keys) (first)) (-> % (vals) (first)) tr/query_ord) tr/pixel_map)
        first_result (first results)]
    (is (= (count results) 10000))
    (is (= (set (keys first_result)) response_set))))

(deftest time-since-change-single-model-test
  (let [result (products/time-since-change (first (:segments tr/first_segments_predictions)) tr/query_ord 100 -100)]
    (is (= (set (keys result)) response_set))))

(deftest time-since-change-chip-level-test
  (let [results (map #(products/time-since-change (-> % (keys) (first)) (-> % (vals) (first)) tr/query_ord) tr/pixel_map)
        greater_thans (filter (fn [i] (> 1000 (:val i))) results)]
    (is (= (count greater_thans) 5887))
    (is (= (count results) 10000))))

(deftest magnitude-of-change-single-model-test
  (let [result (products/magnitude-of-change (first (:segments tr/first_segments_predictions)) tr/query_ord 100 -100)]
    (is (= (set (keys result)) response_set))))

(deftest magnitude-of-change-chip-level-test
  (let [results (map #(products/magnitude-of-change (-> % (keys) (first)) (-> % (vals) (first)) tr/query_ord) tr/pixel_map)
        gt_zero (filter (fn [i] (> (:val i) 0)) results)]
    (is (= (count gt_zero) 387))
    (is (= (count results) 10000))))

(deftest length-of-segment-single-model-test
  (let [result (products/length-of-segment (first (:segments tr/first_segments_predictions)) tr/query_ord 100 -100)]
    (is (= (set (keys result)) response_set))))

(deftest length-of-segment-chip-level-test
  (let [results (map #(products/length-of-segment (-> % (keys) (first)) (-> % (vals) (first)) tr/query_ord) tr/pixel_map)
        gt_zero (filter (fn [i] (> (:val i) 0)) results)]
    (is (= (count gt_zero) 5633))
    (is (= (count results) 10000))))

(deftest curve-fit-single-model-test
  (let [result (products/curve-fit (first (:segments tr/first_segments_predictions)) tr/query_ord 100 -100)]
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
  (let [first_nbr (products/nbr (first tr/first_sorted_segments))
        last_nbr  (products/nbr (last tr/first_sorted_segments))]
    (is (> first_nbr 0.12))
    (is (< first_nbr 0.14))
    (is (> last_nbr  0.33))
    (is (< last_nbr  0.34))))

(deftest get_class_test
  (let [first_class (products/get-class (:prob (first (:predictions tr/first_segments_predictions))))
        last_class  (products/get-class (:prob (last (:predictions tr/first_segments_predictions))))]
    (is (= first_class 7))
    (is (= last_class 5))))

(deftest first_date_of_class_test
  (let [sorted_predictions (util/sort-by-key (:predictions tr/first_segments_predictions) :date)]
    (is (= "1995-07-01" (products/first-date-of-class sorted_predictions 7)))
    (is (= "2012-07-01" (products/first-date-of-class sorted_predictions 5)))
    (is (= nil (products/first-date-of-class sorted_predictions 4)))))

(deftest mean_probabilities_test
  (let [preds [{:prob [0 1 2 3 4 5 6 7 8]} {:prob [7 8 9 5 8 7 6 5 5]}]]
    (is (= [3.5 4.5 5.5 4.0 6.0 6.0 6.0 6.0 6.5]
           (products/mean-probabilities preds)))))

(deftest classify_positive_nbr_test
   (let [model (merge (first tr/first_sorted_segments) {:probabilities tr/grass_to_forest_probs})
         post_forest_query_date (-> "2001-07-01" (util/to-ordinal))
         pre_forest_query_date (-> "1998-07-01" (util/to-ordinal))
         nbrdiff (float 0.06)]
     ; default value 3 is grass, 4 is tree
     (is (= 4 (products/classify model post_forest_query_date 0 nbrdiff)))
     (is (= 3 (products/classify model pre_forest_query_date 0 nbrdiff)))
     (is (= 3 (products/classify model post_forest_query_date 1 nbrdiff)))
     (is (= 4 (products/classify model pre_forest_query_date 1 nbrdiff)))))

(deftest classify_negative_nbr_test
  (let [model (merge (first tr/first_sorted_segments) {:probabilities tr/forest_to_grass_probs})
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
        segment_probabilities (filter (fn [i] (= (util/to-ordinal (:sday i)) sday)) (:predictions tr/first_segments_predictions))
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
  (let [segs_probs tr/first_segments_predictions
        first_seg (merge (first (:segments segs_probs)) {:bday "2001-09-11"})
        last_seg (last (:segments segs_probs))
        modded_segments (merge segs_probs {:segments [first_seg last_seg]})]

    ; query date precedes first segment start date and fill_begin is true
    (is (= (:snow (:lc_map config)) (products/landcover segs_probs (-> "1980-01-01" (util/to-ordinal)) 0)))
    
    ; query date precedes first segment start date
    (is (= (:lc_insuff (:lc_defaults config))
           (products/landcover segs_probs (-> "1980-01-01" (util/to-ordinal)) 0 (merge config {:fill_begin false}))))

    ; query date follows last segment end date and fill_end is true
    (is (= (:water (:lc_map config)) (products/landcover segs_probs (-> "2018-01-01" (util/to-ordinal)) 0)))

    ; query date follows last segment end date
    (is (= (:lc_insuff (:lc_defaults config))
           (products/landcover segs_probs (-> "2017-10-01" (util/to-ordinal)) 0 (merge config {:fill_end false}))))

    ; query date falls between a segments start and end dates
    (is (= (:water (:lc_map config))
           (products/landcover segs_probs (-> "2002-01-01" (util/to-ordinal)) 0)))

    ; query date falls between segments of same landcover classification and fill_samelc config is true
    (is (= (:tree (:lc_map config))
           (products/landcover tr/first_segments_matching_predictions (-> "2001-09-20" (util/to-ordinal)) 0)))

    ; query date falls between one segments break date and the following segments start date and fill_difflc config is true
    (is (= (:water (:lc_map config)) 
           (products/landcover modded_segments (-> "2001-10-03" (util/to-ordinal)) 0)))

    ; query date falls between a segments end date and break date and fill_difflc config is true
    (is (= (:snow (:lc_map config)) 
           (products/landcover modded_segments (-> "2001-09-10" (util/to-ordinal)) 0)))

    ; as a last resort return lc_inbtw configuration value
    (is (= (:lc_inbtw config))
        (products/landcover modded_segments (-> "2001-09-20" (util/to-ordinal)) 0 (merge config {:fill_difflc false})))

    (is (= {:pixelx 1 :pixely 2 :val 75}
           (products/annual-change {:px 1 :py 2} segs_probs (-> "2002-01-01" (util/to-ordinal)))))))


(deftest landcover_confidence_test ; first segment -> sday 1982-12-27 bday 2001-10-04 eday 2001-09-10
                                   ; last segment -> sday 2001-10-04  bday 2017-09-14 eday 2017-09-14
  (let [segs_probs tr/first_segments_predictions
        segs_probs_match_preds tr/first_segments_matching_predictions
        modded_segments (merge segs_probs {:segments [(first (:segments segs_probs))
                                                      (merge (last (:segments segs_probs)) {:chprob 1.0})]})
        ordinal_1995 (-> "1995-07-01" (util/to-ordinal))
        ordinal_2001 (-> "2001-09-20" (util/to-ordinal))
        ordinal_2017 (-> "2017-10-01" (util/to-ordinal))]

    ; query date precedes first segment start date and fill_begin is true
    (is (= (:lcc_back (:lc_defaults config)) 
           (products/confidence segs_probs (-> "1980-01-01" (util/to-ordinal)) 0)))

    ; query date follows last segment end date and change prob is 1
    (is (= (:lcc_afterbr (:lc_defaults config))
           (products/confidence modded_segments ordinal_2017 0)))

     ; query date follows last segment end date and change prob is 0
    (is (= (:lcc_forwards (:lc_defaults config))
           (products/confidence segs_probs ordinal_2017 0)))

    ; query date falls between a segments start date and end date and growth is true
    (is (= (:lcc_growth (:lc_defaults config))
           (products/confidence segs_probs ordinal_1995 0)))

    ; query date falls between a segments start date and end date and decline is true
    (with-redefs [products/nbr (fn [i] -0.66)]
      (is (= (:lcc_decline (:lc_defaults config))
             (products/confidence segs_probs ordinal_1995 0))))

    ; query date falls between a segments start and end date, neither growth nor decline
    (with-redefs [products/nbr (fn [i] 0.01)]
      (is (= 8 (products/confidence segs_probs ordinal_1995 0))))

    ; query date falls between segments of same landcover classification and fill_samelc config is true
    (is (= (:lcc_samelc (:lc_defaults config))
           (products/confidence segs_probs_match_preds ordinal_2001 0)))

    ; query date falls between segments with different landcover classifications
    (is (= (:lcc_difflc (:lc_defaults config))
           (products/confidence segs_probs ordinal_2001 0)))))

(deftest pixel_map_test
  (let [inmap {:pixelxy [1 2] :segments {[1 2] [:a :b :c]} :predictions {[1 2] [:d :e :f]}}]
    (is (= (products/pixel_map inmap) {{:px 1 :py 2} {:segments [:a :b :c] :predictions [:d :e :f]}}))))

(deftest pixel_groups_test
  (let [inmaps [{:px 1 :py 2 :foo "bar"} {:px 1 :py 2 :foo "too"} 
                {:px 3 :py 4 :foo "shizzle"} {:px 3 :py 4 :foo "baz"}]
        value (products/pixel_groups inmaps)] 
    (is (= (count value) 2))
    (is (= (keys value) '([1 2] [3 4])))))

(deftest flatten_product_data_test
  ; lower left pixel value should be the first item in returned collection
  (let [input [{:pixely 3159045, :pixelx -2114775, :val 2} ; lower left
               {:pixely 3159045, :pixelx -2114745, :val 4} ; lower right
               {:pixely 3159075, :pixelx -2114775, :val 6} ; upper left
               {:pixely 3159075, :pixelx -2114745, :val 8} ; upper right
               ]]
    (is (= (products/flatten_product_data input)
           [6 8 2 4]))))

(deftest data_test
  (let [segs  tr/segments_json    ;(:segments tr/first_segments_predictions)  
        preds tr/predictions_json ;(:predictions tr/first_segments_predictions)
        first_seg (:segments tr/first_segments_predictions)
        first_pred (:predictions tr/first_segments_predictions)
        product "time-since-change"
        query_day "2006-07-01"]
    ; requesting a change product with segments and predictions is valid
    (is (= 2915 (nth (products/data segs preds product query_day) 5)))
    ; products/data should throw an exception when its output is too small
    (is (thrown-with-msg? Exception #"Validation Error" (products/data first_seg first_pred product query_day)))
    ; requesting a landcover product with no predictions should throw an exception
    (is (thrown-with-msg? Exception #"Error calculating landcover" (products/data segs [] "primary-landcover" query_day)))
    ; requesting a change product with no predictions is valid
    (is (= 2915 (nth (products/data segs [] product query_day) 5)))))


