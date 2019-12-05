(ns lcmap.gaia.product-specs-test
  (:require [clojure.test :refer :all]
            [clojure.walk :refer [keywordize-keys]]
            [lcmap.gaia.product-specs :as product-specs]
            [lcmap.gaia.test-resources :as tr]))

(deftest segment_check_test
  (let [good_segment (keywordize-keys (first tr/segments_json))
        bad_segment (conj good_segment {:px "not a coordinate"})
        default_segment (conj good_segment {:sday "0001-01-01" :eday "0001-01-01" :bday "0001-01-01"})]
    (is (= good_segment (product-specs/segment_check good_segment)))
    (is (thrown-with-msg? Exception #"Validation Error" (product-specs/segment_check bad_segment)))
    (is (thrown-with-msg? Exception #"Validation Error" (product-specs/segment_check default_segment)))))

(deftest segment_coll_test
  (let [segments (keywordize-keys tr/segments_json)
        validated_segments (product-specs/segment_coll_check segments)]
    (is (= 20166 (count validated_segments)))
    (is (= (first segments) (first validated_segments)))))

(deftest prediction_check_test
  (let [good_prediction (keywordize-keys (first tr/predictions_json))
        bad_prediction_prob (conj good_prediction {:prob "not probabilities"})
        bad_prediction_coord (conj good_prediction {:cx "niner"})]
    (is (= good_prediction (product-specs/prediction_check good_prediction)))
    (is (thrown-with-msg? Exception #"Validation Error" (product-specs/prediction_check bad_prediction_prob)))
    (is (thrown-with-msg? Exception #"Validation Error" (product-specs/prediction_check bad_prediction_coord)))))

(deftest prediction_coll_test
  (let [predictions (keywordize-keys tr/predictions_json)
        validated_predictions (product-specs/prediction_coll_check predictions)]
    (is (= 20166 (count validated_predictions)))
    (is (= (first predictions) (first validated_predictions)))))

