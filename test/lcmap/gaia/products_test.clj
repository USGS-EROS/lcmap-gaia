(ns lcmap.gaia.products-test
  (:require [clojure.test :refer :all]
            [lcmap.gaia.products :as products]
            [lcmap.gaia.file     :as file]
            [lcmap.gaia.util     :as util]))

(def chipdata (file/read-json "resources/y3161805_x-2115585.json"))
(def pixel_segments (util/pixel-groups chipdata))
(def querydate "2006-07-01")

(deftest time-of-change-single-model-test
  (let [first_pixel    (first pixel_segments)
        pixel_models   (last first_pixel)
        result (products/time-of-change (first pixel_models) querydate 100 -100)]
    (is (= (set (keys result))  (set [:pixelx :pixely :toc])))))

(deftest time-of-change-chip-level-test
  (let [results (map #(products/time-of-change (first %) (last %) querydate) pixel_segments)
        first_result (first results)]
    (is (= (count results) 10000))
    (is (= (set (keys first_result)) (set [:pixelx :pixely :toc])))))
