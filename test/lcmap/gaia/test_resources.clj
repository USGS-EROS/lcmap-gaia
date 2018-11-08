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

(def segmentA
  {:thint -12926.188, :bday "2001-10-04", :grmag 112.140785, :blint 4227.913, :grrmse 117.07277, 
   :s2mag 568.2095, :sday "1982-12-27", :grint 14096.05, 
   :grcoef [-0.018689321 -7.062428 52.64909 66.31735 -3.384397 15.727496 17.459877], :px -2114685, 
   :rermse 103.93077, :s1mag 667.05176, :eday "2001-09-10", :nimag 1318.0848, 
   :thcoef [0.01930208 -780.8634 -155.83296 -90.04074 -17.228466 56.151688 -1.1367743], 
   :blcoef [-0.0053713815 37.71041 28.910141 47.774113 -8.368879 36.919247 -4.0322], 
   :s1rmse 104.53059, :nicoef [0.019062052 -242.82335 -177.86938 351.0208 -59.4798 100.0766 56.906433], 
   :s2int 20060.727, :s2rmse 65.64616, :thmag 334.4654, :reint 15451.567, :cx -2115585, :cy 3119805, 
   :s1int 35313.953, :blrmse 131.82863, :remag 257.33282, :thrmse 347.03152, 
   :s2coef [-0.027004562 -22.141205 34.586647 23.296415 -0.0 25.820158 20.58256], :chprob 1.0, 
   :curqa 8, :blmag 78.48317, :niint -10429.385, 
   :s1coef [-0.04701265 -100.99869 2.9391441 68.473 -31.411777 63.906292 18.29616], :py 3118215, 
   :nirmse 233.11906, 
   :recoef [-0.02075027 23.753778 64.33015 39.209644 6.231797 31.255816 0.44192594]})


(def forest_first_prob {:cx -2115585, :cy 3119805, :px -2114685, :py 3118215, :sday "1982-12-27", :eday "2001-09-10", :date "1995-07-01"
                        :prob [0.082641214 0.23013946 1.6417161E-8 0.7288632 3.1727698E-9 7.0137938E-4 0.35765457 6.6413847E-8 6.6413847E-8]})
(def grass_last_prob {:cx -2115585, :cy 3119805, :px -2114685, :py 3118215, :sday "1982-12-27", :eday "2001-09-10", :date "2000-07-01"
                      :prob [6.6413847E-8 6.6413847E-8 0.7137938 0.082641214 0.35765457 3.1727698E-9 1.6417161E-8 0.3288632 0.23013946]})


(def grass_first_prob {:cx -2115585, :cy 3119805, :px -2114685, :py 3118215, :sday "1982-12-27", :eday "2001-09-10", :date "1995-07-01"
                       :prob [6.6413847E-8 6.6413847E-8 0.7137938 0.082641214 0.35765457 3.1727698E-9 1.6417161E-8 0.3288632 0.23013946]})
(def forest_last_prob {:cx -2115585, :cy 3119805, :px -2114685, :py 3118215, :sday "1982-12-27", :eday "2001-09-10", :date "2000-07-01"
                       :prob [0.082641214 0.23013946 1.6417161E-8 0.7288632 3.1727698E-9 7.0137938E-4 0.35765457 6.6413847E-8 6.6413847E-8]})

(def forest_to_grass_probs [forest_first_prob grass_last_prob])
(def grass_to_forest_probs [grass_first_prob forest_last_prob])
(def random_probs [{:cx -2115585, :cy 3119805, :px -2114685, :py 3118215, :sday "1982-12-27", :eday "2001-09-10", :date "1995-07-01"
                    :prob [0.082641214 0.23013946 1.6417161E-8 0.3288632 3.1727698E-9 7.0137938E-4 0.35765457 6.6413847E-8 6.6413847E-8]} 
                   {:cx -2115585, :cy 3119805, :px -2114685, :py 3118215, :sday "1982-12-27", :eday "2001-09-10", :date "2000-07-01"
                    :prob [6.6413847E-8 6.6413847E-8 7.0137938E-4 0.082641214 0.35765457 3.1727698E-9 1.6417161E-8 0.3288632 0.23013946]}] )

(def first_pixelxy    (-> pixel_map (first) (keys) (first)))
(def first_segment    (-> pixel_map (first) (vals) (first) (:segments) (first)))
(def first_prediction (-> pixel_map (first) (vals) (first) (:predictions) (first)))
(def first_probs      (-> first_prediction (:prob)))
(def last_pixelxy     (-> pixel_map (last) (keys) (first)))
(def last_segment     (-> pixel_map (last) (vals) (first) (:segments) (first)))
(def last_prediction  (-> pixel_map (last) (vals) (first) (:predictions) (first)))
(def last_probs       (-> last_prediction (:prob)))
(def grid_data (file/read-json "resources/grid.json"))

