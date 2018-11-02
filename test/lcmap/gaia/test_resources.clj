(ns lcmap.gaia.test-resources
  (:require [lcmap.gaia.file :as file]
            [lcmap.gaia.util :as util]))

(defn ccdc_map
  [inputs]
  (let [pixelx      (first (:pixelxy inputs))
        pixely      (last  (:pixelxy inputs))
        segments    (get (:segments inputs)    (:pixelxy inputs)) 
        predictions (get (:predictions inputs) (:pixelxy inputs))]
    (hash-map {:px pixelx :py pixely} (hash-map :segments segments :predictions predictions))))

(def querydate "2006-07-01")
(def segments    (file/read-json "resources/cx-2115585_cy3119805_segment.json"))
(def predictions (file/read-json "resources/cx-2115585_cy3119805_prediction_with_fake_date.json"))
(def grouped_segments    (group-by (util/variable-juxt ["px" "py"]) segments))
(def grouped_predictions (group-by (util/variable-juxt ["px" "py"]) predictions))
(def first_grouped_predictions (get grouped_predictions [-2114685 3118215]))
(def pixel_map (map #(ccdc_map {:pixelxy % :segments grouped_segments :predictions grouped_predictions}) (keys grouped_segments)))
(def first_pixelxy    (-> pixel_map (first) (keys) (first)))
(def first_segment    (-> pixel_map (first) (vals) (first) (:segments) (first)))
(def first_prediction (-> pixel_map (first) (vals) (first) (:predictions) (first)))
(def first_probs      (-> first_prediction (get "prob")))
(def last_pixelxy     (-> pixel_map (last) (keys) (first)))
(def last_segment     (-> pixel_map (last) (vals) (first) (:segments) (first)))
(def last_prediction  (-> pixel_map (last) (vals) (first) (:predictions) (first)))
(def last_probs       (-> last_prediction (get "prob")))
(def grid_data (file/read-json "resources/grid.json"))

