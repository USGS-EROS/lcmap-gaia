(ns lcmap.gaia.product-specs-test
  (:require [clojure.test :refer :all]
            [clojure.walk :refer [keywordize-keys]]
            [lcmap.gaia.product-specs :as product-specs]
            [lcmap.gaia.test-resources :as tr]))

(deftest segments_check_test
  (let [good_segment (keywordize-keys (first tr/segments_json))
        bad_segment (conj good_segment {:px "not a coordinate"})]
    (is (= good_segment (product-specs/segments_check good_segment)))
    (is (thrown-with-msg? Exception #"Validation Error" (product-specs/segments_check bad_segment)))))

(deftest predictions_check_test
  (let [good_prediction (keywordize-keys (first tr/predictions_json))
        bad_prediction_prob (conj good_prediction {:prob "not probabilities"})
        bad_prediction_coord (conj good_prediction {:cx 123})]
    (is (= good_prediction (product-specs/predictions_check good_prediction)))
    (is (thrown-with-msg? Exception #"Validation Error" (product-specs/predictions_check bad_prediction_prob)))
    (is (thrown-with-msg? Exception #"Validation Error" (product-specs/predictions_check bad_prediction_coord)))
))

