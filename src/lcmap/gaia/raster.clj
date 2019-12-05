(ns lcmap.gaia.raster
  (:require [clojure.tools.logging :as log]
            [clojure.java.io       :as io]
            [clojure.string        :as string]
            [lcmap.gaia.config     :refer [config]]
            [lcmap.gaia.chipmunk   :as chipmunk]
            [lcmap.gaia.file       :as file]
            [lcmap.gaia.gdal       :as gdal]
            [lcmap.gaia.nemo       :as nemo]
            [lcmap.gaia.change-products :as products]
            [lcmap.gaia.storage    :as storage]
            [lcmap.gaia.util       :as util]))

;; a tile is 50 chips x 50 chips
;; each chip is 100 pixels x 100 pixels
;; each pixel is 30m x 30m
(def meters_per_pixel     30)
(def pixels_per_chip_side 100)
(def chips_per_tile_side  50)

;; color definitions for landcover geotiffs 
(def color_nodata    (java.awt.Color.   0   0   0))
(def color_developed (java.awt.Color. 255  50  50))
(def color_cropland  (java.awt.Color. 190 140  90))
(def color_grassland (java.awt.Color. 230 240 210))
(def color_tree      (java.awt.Color.  28  99  48))
(def color_water     (java.awt.Color.   0 112 255))
(def color_wetland   (java.awt.Color. 179 217 255))
(def color_snow      (java.awt.Color. 255 255 255))
(def color_barren    (java.awt.Color. 179 174 163))
(def color_change    (java.awt.Color. 171   0 214))

(def change_vals 
  "Annual Change product values to identify with color_change defined color."
  (filter (fn [i] (not (.contains [11 22 33 44 55 66 77 88] i))) (range 10 89)))

(def change_val_map 
  "Hash-map of change_vals and color_change definition"
  (zipmap change_vals (repeat (count change_vals) color_change)))

(def color_map
  "Associate landcover values with the desired colors"
  (merge 
   (hash-map (:none    (:lc_map config)) color_nodata
             (:develop (:lc_map config)) color_developed
             (:ag      (:lc_map config)) color_cropland
             (:grass   (:lc_map config)) color_grassland
             (:tree    (:lc_map config)) color_tree
             (:water   (:lc_map config)) color_water
             (:wetland (:lc_map config)) color_wetland
             (:snow    (:lc_map config)) color_snow
             (:barren  (:lc_map config)) color_barren)
   change_val_map))

(def product_details
  "Associate product name with its official abbreviation and data type"
  (hash-map "primary-landcover"    {:abbr "LCPRI"   :type gdal/int8}              
            "secondary-landcover"  {:abbr "LCSEC"   :type gdal/int8}            
            "primary-confidence"   {:abbr "LCPCONF" :type gdal/int8}  
            "secondary-confidence" {:abbr "LCSCONF" :type gdal/int8} 
            "annual-change"        {:abbr "LCACHG"  :type gdal/int8} 
            "time-of-change"       {:abbr "SCTIME"  :type gdal/int16}                
            "magnitude-of-change"  {:abbr "SCMAG"   :type gdal/float32}            
            "time-since-change"    {:abbr "SCLAST"  :type gdal/int16}              
            "curve-fit"            {:abbr "SCMQA"   :type gdal/int8}                     
            "length-of-segment"    {:abbr "SCSTAB"  :type gdal/int16}))

(defn get-products
  "Retrieve product details based on type"
  [type]
  (if (= type "cover")
    (filter #(re-matches #"LC(.*)" (:abbr (last %))) product_details)
    (filter #(re-matches #"SC(.*)" (:abbr (last %))) product_details)))

(defn is-landcover
  "Return True if name is a landcover product"
  [name]
  (let [abbr (:abbr (get product_details name))]
    (.contains ["LCPRI" "LCSEC" "LCACHG"] abbr)))

(defn map-details
  "Return attributes of map product"
  [tileid product_info date product]
  (let [grid      (:region config)
        repr_year (first (string/split date #"-"))
        production_date (string/replace (util/todays-date) "-" "") 
        ccd_ver   (:ccd_ver config)
        data_product (first product_info)
        product_abbr (:abbr (last product_info)) 
        elements ["LCMAP" grid tileid repr_year production_date ccd_ver product_abbr]
        name (str (string/join "-" elements) ".tif")
        prefix (storage/get-prefix grid date tileid "raster" product)
        url (storage/get_url storage/dest_bucket (str prefix "/" name))
        type (:type (last product_info))]
    {:name name :prefix prefix :url url :data-type type :data-product data_product}))

(defn calc_offset
  "Returns pixel offset for UL coordinates of chips"
  [tile_x tile_y chip_x chip_y]
  (let [x_diff (- chip_x tile_x)
        y_diff (- tile_y chip_y)
        x_offset (/ x_diff meters_per_pixel)
        y_offset (/ y_diff meters_per_pixel)]
    [(int x_offset) (int y_offset)]))

(defn create_blank_tile_tiff
  "Create empty raster"
  [name ulx uly projection data_type data_product]
  (let [values (repeat (* 5000 5000) 0)
        colortable (if (is-landcover data_product) (gdal/create_colortable color_map) nil)]
      (gdal/create_geotiff name values ulx uly projection data_type 5000 5000 0 0 colortable)))

(defn nlcd_filter
  "Mask, and when appropriate, fill, input data based on nlcd layer"
  ([indata product cx cy filters]
   (let [mask    (:mask filters)
         values  (:values filters)
         fill_fn (fn [input fill] (if (zero? input) fill input))
         masked  (map * indata mask)]
     (if (is-landcover product)
       (map fill_fn masked values)
       masked)))
  ([indata product cx cy]
   (let [filters (chipmunk/nlcd_filters cx cy)]
     (nlcd_filter indata product cx cy filters))))

(defn add_chip_to_tile
  "Add chips worth of data to tile raster"
  [name values tile_x tile_y chip_x chip_y]
  (let [[x_offset y_offset] (calc_offset tile_x tile_y chip_x chip_y)]
    (gdal/update_geotiff name values x_offset y_offset)))

(defn create_raster
  "Create product rasters"
  [{date :date tile :tile tx :tx ty :ty chips :chips product :product :as all}]
  (let [projection   (util/get-projection)
        product_info (get-products product)
        details_fn   #(map-details tile % date product)
        rasters      (into [] (map details_fn product_info))] ; (:name :prefix :url :data-type :data-product)        

    (try
      ; create empty tiffs
      (log/infof "creating empty rasters")
      (doall (map #(create_blank_tile_tiff (:name %) tx ty projection (:data-type %) (:data-product %)) rasters))

      ; step thru the chips coords and retrieve the data for all products
      (log/infof "stepping through chips, adding data to rasters...")
      (doseq [chip chips
              :let [cx (:cx chip)
                    cy (:cy chip)
                    nlcd_data (chipmunk/nlcd_filters cx cy)
                    chip_path (storage/ppath product cx cy tile date)
                    chip_data (util/with-retry (storage/get_json chip_path))
                    chip_vals (fn [i] (nlcd_filter (map #(get % i) chip_data) i cx cy nlcd_data))]]
        
        ; for each raster product, add the appropriate data to the correct raster
        (log/infof "adding cx: %s cy: %s data to rasters" cx cy)
        (doall (map #(add_chip_to_tile (:name %) (chip_vals (:data-product %)) tx ty cx cy) rasters))
        (log/infof "success for cx: %s cy: %s" cx cy))

      (log/infof "done adding chip data to rasters!")
     
      ; push rasters to object store
      (log/infof "pushing tiffs to object storage...")
      (doall (map #(storage/put_tiff % (:name %)) rasters))
      (log/infof "done pushing tiffs")

      ; clean up tiffs
      (log/infof "deleting local tiffs...")
      (doall (map #(io/delete-file (:name %)) rasters))
      (log/infof "done deleting")

      ; return urls
      (doall (map :url rasters)) 

      (catch Exception e
        (let [names (doall (map :name rasters))
              msg (format "problem creating rasters %s: %s" names (.getMessage e))]
          (log/error msg)
          (throw (ex-info msg {:type "data-generation-error" :message msg})))))))
