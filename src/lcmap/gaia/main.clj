(ns lcmap.gaia.main
  (:gen-class)
  (:require [clojure.tools.logging :as log]
            [lcmap.gaia.config     :as config]
            [lcmap.gaia.file       :as file]
            [lcmap.gaia.products   :as products]
            [lcmap.gaia.util       :as util]
            [mount.core            :as mount]
            [lcmap.gaia.gdal       :as gdal]))



;; ccdc results are the number of change segments
;; detected over a 100x100 pixel area (chip)
;; There can be multiple changes over time 
;; detected.
;; each pixel will have at least 1 entry
;; if no changes are detected, a result is 
;; returned with an sday and eday == 0


(defmulti gen-product
  (fn [infile product queryday] (keyword product)))

(defmethod gen-product :default
  [infile product queryday]
  (log/infof "%s is not a valid product" product)
  nil)

(defmethod gen-product :time-of-change
  [infile product queryday]
  (let [input (file/read-json infile)
        pixel_segments (util/pixel-groups input)
        pixel_array (map #(products/time-of-change (first %) (last %) queryday) pixel_segments)
        output_name (products/product-name (first input) product "gtif")]
    ;(gdal/geotiff pixel_array output_name)
    output_name))

(defn -main
  ([]
   (println "Welcome to lcmap-gaia\nHere are your options:")
   (System/exit 0))
  ([infile product queryday]
   ;; Only mount states defined in required namespaces are started.
   (mount/start)
   
   (log/infof "output for %s %s %s" infile product)
   (System/exit 0)))

