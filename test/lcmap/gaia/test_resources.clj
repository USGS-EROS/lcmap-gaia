(ns lcmap.gaia.test-resources
  (:require [lcmap.gaia.file :as file]
            [lcmap.gaia.util :as util]
            [clojure.walk :refer [keywordize-keys]]))

(defn ccdc_map
  [inputs]
  (let [pixelx      (first (:pixelxy inputs))
        pixely      (last  (:pixelxy inputs))
        segments    (get (:segments inputs)    (:pixelxy inputs)) 
        predictions (get (:predictions inputs) (:pixelxy inputs))]
    (hash-map {:px pixelx :py pixely} (hash-map :segments segments :predictions predictions))))

(def querydate "2006-07-01")
(def query_ord (util/to-ordinal querydate))
(def segments    (file/read-json "resources/cx-2115585_cy3119805_segment.json"))
(def predictions (file/read-json "resources/cx-2115585_cy3119805_prediction_with_fake_date.json"))
(def grouped_segments    (-> ["px" "py"] (util/variable-juxt) (group-by segments) (keywordize-keys)))
(def grouped_predictions (-> ["px" "py"] (util/variable-juxt) (group-by predictions) (keywordize-keys)))
(def first_grouped_predictions (get grouped_predictions [-2114685 3118215]))
(def pixel_map (map #(ccdc_map {:pixelxy % :segments grouped_segments :predictions grouped_predictions}) (keys grouped_segments)))
(def first_pixelxy (-> (first pixel_map) (keys) (first)))
(def first_segments_predictions (-> (first pixel_map) (vals) (first)))
(def first_sorted_segments (util/sort-by-key (:segments first_segments_predictions) :sday))
(def first_probabilities (:predictions first_segments_predictions))



(def first_pixelxy    (-> pixel_map (first) (keys) (first)))
(def first_segment    (-> pixel_map (first) (vals) (first) (:segments) (first)))
(def first_prediction (-> pixel_map (first) (vals) (first) (:predictions) (first)))
(def first_probs      (-> first_prediction (:prob)))
(def last_pixelxy     (-> pixel_map (last) (keys) (first)))
(def last_segment     (-> pixel_map (last) (vals) (first) (:segments) (first)))
(def last_prediction  (-> pixel_map (last) (vals) (first) (:predictions) (first)))
(def last_probs       (-> last_prediction (:prob)))
(def grid_data (file/read-json "resources/grid.json"))

