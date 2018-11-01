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
(def pixel_map (map #(ccdc_map {:pixelxy % :segments grouped_segments :predictions grouped_predictions}) (keys grouped_segments)))
(def grid_data (file/read-json "resources/grid.json"))
