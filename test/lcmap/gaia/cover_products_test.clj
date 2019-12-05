(ns lcmap.gaia.cover-products-test
  (:require [clojure.test :refer :all]
            [clojure.string      :as string]
            [clojure.walk :refer [stringify-keys]]
            [clojure.math.combinatorics :as combo]
            [lcmap.gaia.cover-products :as cover-products]
            [lcmap.gaia.file     :as file]
            [lcmap.gaia.util     :as util]
            [lcmap.gaia.config   :refer [config]]
            [lcmap.gaia.storage  :as storage]
            [lcmap.gaia.test-resources :as tr]))


(deftest falls_between_eday_sday-coll-test
  (let [map_a {:follows_eday true :precedes_sday false}
        map_b {:precedes_sday true :follows_eday false}
        expected [map_a map_b]]
    (is (= (cover-products/falls-between-eday-sday map_a map_b) expected))))

(deftest falls_between_eday_sday-map-test
  (let [map_a {:follows_eday false :precedes_sday true}
        map_b {:precedes_sday true :follows_eday false}]
    (is (= (cover-products/falls-between-eday-sday map_a map_b) map_b))))

(deftest falls_between_eday_sday-nonmap-test
  (let [map_a {:follows_eday true :precedes_sday false}
        map_b {:precedes_sday true :follows_eday false}
        map_c {:precedes_sday true :follows_eday false}
        expected [map_a map_b]]
    (is (= (cover-products/falls-between-eday-sday [map_a map_b] map_c) expected))))

(deftest falls_between_bday_sday-coll-test
  (let [map_a {:follows_bday true :precedes_sday false}
        map_b {:precedes_sday true :follows_bday false}
        expected [map_a map_b]]
    (is (= (cover-products/falls-between-bday-sday map_a map_b) expected))))

(deftest nbr_test
  (let [first_seg (first tr/first_sorted_segments)
        first_sday (util/to-ordinal (get first_seg "sday"))
        first_eday (util/to-ordinal (get first_seg "eday"))
        last_seg (last tr/first_sorted_segments)
        last_sday (util/to-ordinal (get last_seg "sday"))
        last_eday (util/to-ordinal (get last_seg "eday"))
        first_nbr (cover-products/normalized-burn-ratio first_seg first_sday first_eday)
        last_nbr  (cover-products/normalized-burn-ratio last_seg last_sday last_eday)]
    (is (> first_nbr 0.12))
    (is (< first_nbr 0.14))
    (is (> last_nbr  0.33))
    (is (< last_nbr  0.34))))

(deftest get_class_test
  (let [first_class (cover-products/get-class (get (first (:predictions tr/first_segments_predictions)) "prob"))
        last_class  (cover-products/get-class (get (last (:predictions tr/first_segments_predictions)) "prob"))]
    (is (= first_class 6))
    (is (= last_class 4))))

(deftest first_date_of_class_test
  (let [sorted_predictions (util/sort-by-key (:predictions tr/first_segments_predictions) "pday")]
    (is (= "1995-07-01" (cover-products/first-date-of-class sorted_predictions 6)))
    (is (= "2012-07-01" (cover-products/first-date-of-class sorted_predictions 4)))
    (is (= nil (cover-products/first-date-of-class sorted_predictions 3)))))

(deftest mean_probabilities_test
  (let [preds [{"prob" [0 1 2 3 4 5 6 7 8]} {"prob" [7 8 9 5 8 7 6 5 5]}]]
    (is (= [3.5 4.5 5.5 4.0 6.0 6.0 6.0 6.0 6.5]
           (cover-products/mean-probabilities preds)))))

(deftest classify_positive_nbr_test
   (let [predictions tr/grass_to_forest_probs
         segments tr/first_sorted_segments
         post_forest_query_date (-> "2001-07-01" (util/to-ordinal))
         pre_forest_query_date (-> "1998-07-01" (util/to-ordinal))
         nbrdiff (float 0.06)]
     ;; default value 3 is grass, 4 is tree
     (is (= 4 (cover-products/classify predictions post_forest_query_date 0 nbrdiff)))
     (is (= 3 (cover-products/classify predictions pre_forest_query_date 0 nbrdiff)))
     (is (= 3 (cover-products/classify predictions post_forest_query_date 1 nbrdiff)))
     (is (= 4 (cover-products/classify predictions pre_forest_query_date 1 nbrdiff)))
     ))

(deftest classify_negative_nbr_test
  (let [predictions tr/forest_to_grass_probs  ;(map cover-products/convert_prediction_dates tr/forest_to_grass_probs)
        segments tr/first_sorted_segments
        post_grass_query_date (-> "2001-07-01" (util/to-ordinal))
        pre_grass_query_date (-> "1998-07-01" (util/to-ordinal))
        nbrdiff (float -0.06)]
    ;; default value 3 is grass, 4 is tree
    (is (= 3 (cover-products/classify predictions post_grass_query_date 0 nbrdiff)))
    (is (= 4 (cover-products/classify predictions pre_grass_query_date 0 nbrdiff)))
    (is (= 4 (cover-products/classify predictions post_grass_query_date 1 nbrdiff)))
    (is (= 3 (cover-products/classify predictions pre_grass_query_date 1 nbrdiff)))))

(deftest classify_else_test
  (let [first_segment (first tr/first_sorted_segments)
        sday (util/to-ordinal (get first_segment "sday"))
        eday (util/to-ordinal (get first_segment "eday"))
        nbrdiff (cover-products/normalized-burn-ratio first_segment sday eday) ; 0.1204...
        probs (:predictions tr/first_segments_predictions) ; (map cover-products/convert_prediction_dates (:predictions tr/first_segments_predictions))
        sorted_probabilities (util/sort-by-key probs "pday")
        ]
    (is (= 3 (cover-products/classify sorted_probabilities tr/query_ord 0 nbrdiff)))))

(deftest characterize_segment_test
  (with-redefs [cover-products/normalized-burn-ratio (fn [i x z] 66)
                cover-products/classify (fn [a b c d] 99)]
    (let [segment   {"sday" "1990-04-27" "eday" "2000-06-11" "bday" "2000-06-11" "chprob" 1.0}
          query_day (-> "1998-07-01" (util/to-ordinal))
          probabilities [{"sday" "1990-04-27" "pday" "1998-07-01"} 
                         {"sday" "1990-04-27" "pday" "2001-07-01"}]
          characterized (cover-products/characterize-segment segment query_day probabilities)]
      (is (= characterized {:intersects true
                            :precedes_sday false
                            :follows_eday false
                            :follows_bday false
                            :btw_eday_bday false
                            :sday (util/to-ordinal "1990-04-27") 
                            :eday (util/to-ordinal "2000-06-11")  
                            :bday (util/to-ordinal "2000-06-11") 
                            :growth false
                            :decline false
                            :chprob 1.0
                            :probabilities [{"sday" "1990-04-27", "pday" "1998-07-01"} {"sday" "1990-04-27" "pday" "2001-07-01"}]
                            :primary_class 99
                            :secondary_class 99})))))

(deftest landcover_test
  (let [pixel_input tr/pixel_input
        ordinal_date (-> "2001-07-01" util/to-ordinal)
        pixel_dates (combo/cartesian-product [ordinal_date] (keys pixel_input))
        characterized_pixels (map #(cover-products/characterize-inputs (last %) (get pixel_input (last %)) (first %)) pixel_dates)
        first_pixel (first characterized_pixels) ; dates 730301 and 730666 (2001-07-01)
        first_pixel_date (first (keys (:current first_pixel)))
        all_segments (first_pixel_date (:current first_pixel))
        first_segment (first (first_pixel_date (:current first_pixel)))
        first_predictions (:probabilities first_segment)]

    ;; query date precedes first segment start date and fill_begin is true
    ;; query date first_pixel 730666
    ;; first_segment sday 724514
    ;; (:fill_begin config) is true by default
    ;; make query_date less than 724514
    (is (= (:ag (:lc_map config))
           (:primary (cover-products/landcover all_segments 724114))
           ;(cover-products/landcover all_segments 724114 0)
))

    ;; query date follows last segment end date and fill_end is true
    ;; first_segment eday 736594
    ;; (:fill_end config) defaults true
    (is (= (:ag (:lc_map config))
           (:primary (cover-products/landcover all_segments 736894))
           ;(cover-products/landcover all_segments 736894 0)
           ))

    ;; query date falls between a segments start and end dates
    ;; sday 724514, eday 736594
    ;; first_segment intersects is true
    (is (= (:ag (:lc_map config))
           (:primary (cover-products/landcover all_segments 727514))
           ;(cover-products/landcover all_segments 727514 0)
           ))

    ;; query date falls between segments of same landcover classification and fill_samelc config is true
    ;; fill_samelc defaults to true
    (let [segment1 (merge first_segment {:eday 725000 :intersects false :follows_eday true})
          segment2 (merge first_segment {:sday 727000 :intersects false :precedes_sday true})
          input [segment1 segment2] ]
      (is (= (:ag (:lc_map config))
             (:primary (cover-products/landcover input 726000))
             ;(cover-products/landcover input 726000 0)
             )))

    ;; query date falls between one segments break date and the following segments start date and fill_difflc config is true
    ;; fill_difflc defaults true
    (let [segment1 (merge first_segment {:eday 725000 :bday 72500 :intersects false :follows_bday true})
          segment2 (merge first_segment {:sday 727000 :intersects false :precedes_sday true})
          input [segment1 segment2]]
      (is (= (:ag (:lc_map config))
             (:primary (cover-products/landcover input 726000 (merge config {:fill_samelc false})))
             ;(cover-products/landcover input 726000 0 (merge config {:fill_samelc false}))
             )))

    ;; query date falls between a segments end date and break date and fill_difflc config is true
    ;; first_segment: bday 736594, sday 724514, eday 736594
    ;; default date 730666 (2001-07-01)
    ;; fill_difflc config defaults true
    (let [segment1 (merge first_segment {:bday 736894 :intersects false :btw_eday_bday true})
          segment2 (merge first_segment {:sday 736994 :eday 737001 :bday 737001 :intersects false :precedes_sday true})
          input [segment1 segment2] ; query date falls btw eday and bday of 1st segment
          mod_cfg (merge config {:fill_end false :fill_samelc false})]
      (is (= (:ag (:lc_map config))
             (:primary (cover-products/landcover input 736694 mod_cfg))
            ; (cover-products/landcover input 736694 0 mod_cfg)
             ))

      ;; pixel is unclassifiable, throw and exception!  clojure.lang.ExceptionInfo: problem calculating landcover with pixel
      (is (thrown-with-msg? Exception #"problem calculating landcover"
                            (:primary (cover-products/landcover input 736694 (merge mod_cfg {:fill_difflc false})))
                            ;(cover-products/landcover input 736694 0 (merge mod_cfg {:fill_difflc false}))
                            )))

    ))

(deftest landcover_confidence_test
  (let [pixel_input tr/pixel_input
        ordinal_date (-> "2001-07-01" util/to-ordinal)
        pixel_dates (combo/cartesian-product [ordinal_date] (keys pixel_input))
        characterized_pixels (map #(cover-products/characterize-inputs (last %) (get pixel_input (last %)) (first %)) pixel_dates)
        first_pixel (first characterized_pixels) ; dates 730301 and 730666 (2001-07-01)
        first_pixel_date (first (keys (:current first_pixel)))
        all_segments (first_pixel_date (:current first_pixel))
        first_segment (first (first_pixel_date (:current first_pixel)))
        first_predictions (:probabilities first_segment)]

    ;; query date precedes first segment start date and fill_begin is true
    (is (= (:lcc_back (:lc_defaults config))
           (:primary (cover-products/confidence all_segments 724114))
           ;(cover-products/confidence all_segments 724114 0)
           ))

    ;; query date follows last segment end date and change prob is 1
    (is (= (:lcc_afterbr (:lc_defaults config))
           (:primary (cover-products/confidence [(merge first_segment {:chprob 1.0})] 736894))
           ;(cover-products/confidence [(merge first_segment {:chprob 1.0})] 736894 0)
           ))

    ;; query date follows last segment end date and change prob is 0
    (is (= (:lcc_forwards (:lc_defaults config))
           (:primary (cover-products/confidence all_segments 736894))
           ;(cover-products/confidence all_segments 736894 0)
           ))

    ;; query date falls between a segments start date and end date and growth is true
    ; first_pixel growth is true
    (is (= 82 ;(:lcc_growth (:lc_defaults config))
           (:primary (cover-products/confidence all_segments 727514))
           ;(cover-products/confidence all_segments 727514 0)
           ))

    ;; query date falls between a segments start date and end date and decline is true
    (is (= (:lcc_decline (:lc_defaults config))
           (:primary (cover-products/confidence [(merge first_segment {:growth false :decline true})] 727514))
           ;(cover-products/confidence [(merge first_segment {:growth false :decline true})] 727514 0)
           ))

    ;; query date falls between a segments start and end date, neither growth nor decline
    (is (= 82 
           (:primary (cover-products/confidence [(merge first_segment {:growth false})] 727514))
           ;(cover-products/confidence [(merge first_segment {:growth false})] 727514 0)
           ))

    ;; query date falls between segments of same landcover classification and fill_samelc config is true
    (let [segment1 (merge first_segment {:eday 725000 :intersects false :follows_eday true})
          segment2 (merge first_segment {:sday 727000 :intersects false :precedes_sday true})
          input [segment1 segment2]]
      (is (= (:lcc_samelc (:lc_defaults config))
             (:primary (cover-products/confidence input 726000))
             ;(cover-products/confidence input 726000 0)
             )))

    ;; query date falls between segments with different landcover classifications
    (let [segment1 (merge first_segment {:eday 725000 :intersects false :follows_eday true})
          segment2 (merge first_segment {:sday 727000 :intersects false :precedes_sday true :primary_class 3})
          input [segment1 segment2]]
      (is (= (:lcc_difflc (:lc_defaults config))
             (:primary (cover-products/confidence input 726000))
             ;(cover-products/confidence input 726000 0)
             )))
))

;; annual change specific test
(deftest annual-change
  (let [cover (tr/cover_data_1086765_1975335)]
    ;(is (= 41 (:annual-change (:values (first cover)))))
    ))

(deftest generate-test
  (with-redefs [storage/segments-sorted (fn [a b c] (util/sort-by-key tr/segments_json "sday"))
                storage/predictions (fn [a b] tr/predictions_json)
                storage/put_json (fn [a b] true)]
    (let [params {:dates ["2007-07-01" "2008-07-01"] :cx 111111 :cy 222222 :tile "123456"}
          result (cover-products/generate params)]
      (is (= result {:products "cover" :cx 111111 :cy 222222 :dates ["2007-07-01" "2008-07-01"]})))))


