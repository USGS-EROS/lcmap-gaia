(ns lcmap.gaia.main
  (:gen-class)
  (:require [clojure.tools.logging :as log]
            [lcmap.gaia.config     :refer [config]]
            [lcmap.gaia.file       :as file]
            [lcmap.gaia.products   :as products]
            [lcmap.gaia.util       :as util]
            [mount.core            :as mount]
            [lcmap.gaia.gdal       :as gdal]))

(defn flatten-vals
  "Flatten the values for a collection of hash-maps"
  [coll mapkey]
  (let [coll_vals (map (fn [i] (vals i)) coll)
        vals_flat (flatten coll_vals)]
    (map mapkey vals_flat)))

(defn product-data
  "Returns a flat list of product values from JSON of a chips worth of CCDC results"
  [injson product queryday]
  (let [product_fn (resolve (symbol (str "products/" product)))
        ; group segments by pixel coordinates
        pixel_segments (util/coll-groups injson ["pixelx" "pixely"])
        ; map the products function across the pixel segments. Returns a flat
        ; collection, one hash map per pixel coordinate pair.
        pixel_array (map #(product_fn (first %) (last %) queryday) pixel_segments)
        ; group product coll by row 
        ; [{:pixely 3159045} [{:pixely 3159045, :pixelx -2114775, :val 6290},...] ...]
        row_groups (util/coll-groups pixel_array [:pixely]) 
        ; sort row group values by pixelx descending 
        sort-pixelx-fn (fn [i] (hash-map (:pixely (first i)) (sort-by :pixelx (last i))))
        sorted-x-vals (map sort-pixelx-fn row_groups)
        ; sort the rows by the pixely key ascending
        sorted-y-rows (sort-by (fn [i] (first (keys i))) sorted-x-vals)]
    ; finally, flatten to a one dimensional list
    (flatten-vals sorted-y-rows :val)))

;; ccdc results are the number of change segments
;; detected over a 100x100 pixel area (chip)
;; There can be multiple changes over time 
;; detected per pixel.
;; if no changes are detected, a result is 
;; returned with an sday and eday == 0

(defn generate-product
  [infile product queryday]
  (let [input (file/read-json infile)
        pixel_vals_flat (product-data input product queryday)
        output_name (products/product-name (first input) product "tif")
        first_model (first (last pixel_segments))
        chipx (get first_model "chipx")
        chipy (get first_model "chipy")
        proj_wkt (util/get-projection)]
    (gdal/geotiff_from_pixel_array pixel_vals_flat output_name chipx chipy proj_wkt)
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

