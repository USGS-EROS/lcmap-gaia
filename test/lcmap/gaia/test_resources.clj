(ns lcmap.gaia.test-resources
  (:require [lcmap.gaia.file :as file]
            [lcmap.gaia.util :as util]
            [clojure.walk :refer [keywordize-keys]]))

(defn ccdc_map
  [inputs]
  (let [segments    (get (:segments inputs)    (:pixelxy inputs)) 
        predictions (get (:predictions inputs) (:pixelxy inputs))]
    (hash-map {:px (first (:pixelxy inputs)) :py (last (:pixelxy inputs))} 
              {:segments segments :predictions predictions})))

(def query_ord (util/to-ordinal "2006-07-01"))

(def segments_json (file/read-json "resources/cx-2115585_cy3119805_segment.json"))
(def predictions_json (file/read-json "resources/cx-2115585_cy3119805_prediction_with_fake_date.json"))

(def grouped_segments    (-> ["px" "py"] (util/variable-juxt) (group-by segments_json) (keywordize-keys)))
(def grouped_predictions (-> ["px" "py"] (util/variable-juxt) (group-by predictions_json) (keywordize-keys)))
(def pixel_map (map #(ccdc_map {:pixelxy % :segments grouped_segments :predictions grouped_predictions}) (keys grouped_segments)))

(def first_segments_predictions (-> (first pixel_map) (vals) (first)))
(def first_sorted_segments (util/sort-by-key (:segments first_segments_predictions) :sday))

(def first_prob  {:cx -2115585, :cy 3119805, :px -2114685, :py 3118215, :sday "1982-12-27", :eday "2001-09-10", :date "1995-07-01"})
(def last_prob   {:cx -2115585, :cy 3119805, :px -2114685, :py 3118215, :sday "1982-12-27", :eday "2001-09-10", :date "2000-07-01"})
(def forest_prob {:prob [0.082641214 0.23013946 1.6417161E-8 0.7288632 3.1727698E-9 7.0137938E-4 0.35765457 6.6413847E-8 6.6413847E-8]})
(def grass_prob  {:prob [6.6413847E-8 6.6413847E-8 0.7137938 0.082641214 0.35765457 3.1727698E-9 1.6417161E-8 0.3288632 0.23013946]})
(def forest_to_grass_probs [(merge first_prob forest_prob) (merge last_prob grass_prob)])
(def grass_to_forest_probs [(merge first_prob grass_prob) (merge last_prob forest_prob)])
(def first_segments_matching_predictions (merge first_segments_predictions {:predictions [(merge first_prob forest_prob) 
                                                                                          (merge last_prob {:sday "2001-10-04", :eday "2017-09-14", :date "2017-09-14"} forest_prob)]}))



