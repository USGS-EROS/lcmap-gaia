(ns lcmap.gaia.main
  (:gen-class)
  (:require [clojure.tools.logging :as log]
            [cheshire.core         :as json]
            [lcmap.gaia.config     :refer [config]]
            [lcmap.gaia.file       :as file]
            [lcmap.gaia.products   :as products]
            [lcmap.gaia.server     :as server]))

;; ccdc results are the number of change segments
;; detected over a 100x100 pixel area (chip)
;; There can be multiple changes over time 
;; detected per pixel.
;; if no changes are detected, a result is 
;; returned with an sday and eday == 0

(defn -main
  ([]
   (server/run-server))
  ([segments_file  predictions_file product queryday]
   (let [segments_input (file/read-json segments_file)
         predictions_input (file/read-json predictions_file)
         product_data (products/data segments_input predictions_input product queryday)
         chipx (get (first segments_input) "cx")
         chipy (get (first segments_input) "cy")]
     (println (json/generate-string {:x chipx :y chipy :values product_data})))))

