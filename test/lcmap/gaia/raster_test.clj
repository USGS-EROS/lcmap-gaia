(ns lcmap.gaia.raster-test
  (:require [clojure.test        :refer :all]
            [clojure.java.io     :as io]
            [lcmap.gaia.gdal     :as gdal]
            [lcmap.gaia.chipmunk :as chipmunk]
            [lcmap.gaia.products :as products]
            [lcmap.gaia.storage  :as storage]
            [lcmap.gaia.util     :as util]
            [lcmap.gaia.raster   :as raster]))

(deftest test-calc_offset
  (let [result (raster/calc_offset 1500 2400 1560 2370)]
    (is (= result [2 1]))))

(deftest test-create_blank_tile_tiff
  (with-redefs [gdal/create_geotiff (fn [a b c d e f g h i] {:name a :values b :ulx c :uly d :projection e :xdim f :ydim g :xoff h :yoff i})]
    (let [result (raster/create_blank_tile_tiff "foo.tif" 111 222 "wkt_proj")]
      (is (= (sort '(:name :values :ulx :uly :projection :xdim :ydim :xoff :yoff)) (sort (keys result))))
      (is (= (count (:values result)) 25000000))
      (is (= 0 (first (:values result)))))))

(deftest test-add_chip_to_tile
  (with-redefs [gdal/update_geotiff (fn [a b c d] {:name a :values b :xoff c :yoff d})]
    (let [result (raster/add_chip_to_tile "foo.tif" [1 2 3] 300 600 330 570)]
      (is (= result {:name "foo.tif" :values [1 2 3] :xoff 1 :yoff 1})))))

(deftest test-nlcd_filter
  (with-redefs [chipmunk/nlcd_filters (fn [a b] {:mask [0 1 0 1] :values [1 2 3 4]})]
    (let [result (raster/nlcd_filter [5 6 7 0] "time-since-change" 111 222)]
      (is (= result [0 6 0 0])))))

(deftest test-create_geotiff
  (with-redefs [util/get-projection           (fn [] "wkt-proj")
                products/map-path             (fn [a b c] (str a b c))
                raster/create_blank_tile_tiff (fn [a b c d] true)
                products/ppath                (fn [a b c d e] (str a b c d e))
                storage/get_json              (fn [i] {"values" [1 2 3 4]})
                raster/nlcd_filter            (fn [a b c d] [6 7 8 9])
                raster/add_chip_to_tile       (fn [a b c d e f] true)
                storage/put_tiff              (fn [a b] true)
                io/delete-file                (fn [i] true)]
    (let [input {:date "2007-01-01" :tile "001002" :tilex 111111 :tiley 222222 
                 :chips [{:cx 1 :cy 2 :value 3} {:cx 1 :cy 3 :value 4}] 
                 :product "time-since-change"}
          result (raster/create_geotiff input)]
      (is (= result "001002time-since-change2007-01-01")))))

