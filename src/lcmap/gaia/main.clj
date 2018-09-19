(ns lcmap.gaia.main
  (:gen-class)
  (:require [clojure.tools.logging :as log]
            [mount.core            :as mount]
            [lcmap.gaia.config     :refer [config]]
            [lcmap.gaia.products   :as products]
            [lcmap.gaia.gdal       :as gdal]
            [lcmap.gaia.server     :as server]))

;; ccdc results are the number of change segments
;; detected over a 100x100 pixel area (chip)
;; There can be multiple changes over time 
;; detected per pixel.
;; if no changes are detected, a result is 
;; returned with an sday and eday == 0

(defn -main
  ([]
   (mount/start)
   (server/run-server))
  ([infile product queryday]
   ;; Only mount states defined in required namespaces are started.
   (mount/start)
   (products/generate-product infile product queryday)
   (System/exit 0)))

