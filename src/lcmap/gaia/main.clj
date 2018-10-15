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
  ([infile product queryday]
   (let [input (file/read-json infile)
         product_data (products/data input product queryday)
         chipx (get (first input) "chipx")
         chipy (get (first input) "chipy")]
     (println (json/generate-string {:x chipx :y chipy :values product_data})))))

