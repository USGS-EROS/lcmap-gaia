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
;; detected.
;; each pixel will have at least 1 entry
;; if no changes are detected, a result is 
;; returned with an sday and eday == 0

(defn generate-product
  [infile product queryday]
  (let [input (file/read-json infile)
        product_fn (resolve (symbol (str "products/" product)))
        pixel_segments (util/pixel-groups input)
        pixel_array (map #(product_fn (first %) (last %) queryday) pixel_segments)
        output_name (products/product-name (first input) product "gtif")
        first_model (first (last pixel_segments))
        chipx (get first_model "chipx")
        chipy (get first_model "chipy")
        proj_wkt (util/get-projection)]
    ;(gdal/geotiff_from_pixel_array pixel_array output_name chipx chipy proj_wkt)
    (prn output_name)
    output_name))

(defn -main
  ([]
   (println "Welcome to lcmap-gaia\nHere are your options:")
   (System/exit 0))
  ([infile product queryday]
   ;; Only mount states defined in required namespaces are started.
   (mount/start)
   (generate-product infile product queryday)
   (System/exit 0)))

