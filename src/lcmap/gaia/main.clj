(ns lcmap.gaia.main
  (:gen-class)
  (:require [clojure.tools.logging :as log]
            [lcmap.gaia.config     :refer [config]]
            [lcmap.gaia.file       :as file]
            [lcmap.gaia.products   :as products]
            [lcmap.gaia.util       :as util]
            [mount.core            :as mount]
            [lcmap.gaia.gdal       :as gdal]))

;; ccdc results are the number of change segments
;; detected over a 100x100 pixel area (chip)
;; There can be multiple changes over time 
;; detected per pixel.
;; if no changes are detected, a result is 
;; returned with an sday and eday == 0

(defn generate-product
  [infile product queryday]
  (let [input (file/read-json infile)
        product_fn (resolve (symbol (str "lcmap.gaia.products/" product)))
        product_values (products/data input product_fn queryday)
        output_name (products/product-name (first input) product "tif")
        chipx (get (first input) "chipx")
        chipy (get (first input) "chipy")
        proj_wkt (util/get-projection)]
    (gdal/geotiff_from_pixel_array product_values output_name chipx chipy proj_wkt)
    (prn output_name)
    true))

(defn -main
  ([]
   (println "Welcome to lcmap-gaia\nHere are your options:")
   (System/exit 0))
  ([infile product queryday]
   ;; Only mount states defined in required namespaces are started.
   (mount/start)
   (generate-product infile product queryday)
   (System/exit 0)))

