(ns lcmap.gaia.raster
  (:require [clojure.tools.logging :as log]
            [clojure.java.io :as io]
            [lcmap.gaia.gdal :as ggdal]
            [lcmap.gaia.products :as products]
            [lcmap.gaia.nemo :as nemo]
            [lcmap.gaia.file :as file]
            [lcmap.gaia.gdal :as gdal]
            [lcmap.gaia.util :as util]
            [lcmap.gaia.storage :as storage]
            [lcmap.gaia.chipmunk :as chipmunk]))

;; a tile is 50 chips x 50 chips
;; each chip is 100 pixels x 100 pixels
;; each pixel is 30m x 30m
(def meters_per_pixel     30)
(def pixels_per_chip_side 100)
(def chips_per_tile_side  50)

(defn calc_offset
  "Returns pixel offset for UL coordinates of chips"
  [tile_x tile_y chip_x chip_y]
  (let [x_diff (- chip_x tile_x)
        y_diff (- tile_y chip_y)
        x_offset (/ x_diff meters_per_pixel)
        y_offset (/ y_diff meters_per_pixel)]
    [(int x_offset) (int y_offset)]))

(defn create_blank_tile_tiff
  [name ulx uly projection]
  (let [values (repeat (* 5000 5000) 0)]
    (gdal/create_geotiff name values ulx uly projection 5000 5000 0 0)))

(defn add_chip_to_tile
  [name values tile_x tile_y chip_x chip_y]
  (let [[x_offset y_offset] (calc_offset tile_x tile_y chip_x chip_y)]
    (gdal/update_geotiff name values x_offset y_offset)))

(defn nlcd_filter
  [indata product cx cy]
  (let [filters (chipmunk/nlcd_filters cx cy)
        mask (:mask filters)
        values (:values filters)]
    (map * indata mask)))

(defn create_geotiff
  [{date :date tile :tile tilex :tilex tiley :tiley chips :chips product :product :as all}]
  (let [projection (util/get-projection)
        map_path (products/map-path tile product date)]
    (create_blank_tile_tiff (:name map_path) tilex tiley projection)
    (doseq [chip chips
            :let [cx (:cx chip)
                  cy (:cy chip)
                  chip_path (products/ppath product cx cy tile date)
                  chip_data (storage/get_json chip_path)
                  chip_vals (nlcd_filter (get chip_data "values") product cx cy)]]
      (log/debugf "adding %s to tile: %s" (:name chip_path) tile)
      (if chip_data
        (add_chip_to_tile (:name map_path) chip_vals tilex tiley cx cy)
        (log/debugf "no data to add to tile %s at cx: %s | cy: %s" tile cx cy)))
    (log/infof "pushing tiff to object storage: %s" map_path)
    (storage/put_tiff map_path (:name map_path))
    (log/infof "deleting local tiff: %s" (:name map_path))
    (io/delete-file (:name map_path))
    map_path))


