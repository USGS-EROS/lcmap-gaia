(ns lcmap.gaia.main
  (:gen-class)
  (:require [clojure.tools.logging :as log]
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
   (products/data (file/read-json infile) product queryday)
   (System/exit 0)))

