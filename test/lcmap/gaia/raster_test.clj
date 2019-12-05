(ns lcmap.gaia.raster-test
  (:require [clojure.test        :refer :all]
            [clojure.java.io     :as io]
            [clojure.string      :as string]
            [clojure.tools.logging :as log]
            [lcmap.gaia.gdal     :as gdal]
            [lcmap.gaia.chipmunk :as chipmunk]
            [lcmap.gaia.storage  :as storage]
            [lcmap.gaia.util     :as util]
            [lcmap.gaia.raster   :as raster]
            [lcmap.gaia.config :refer [config]]))

(deftest test-calc_offset
  (let [result (raster/calc_offset 1500 2400 1560 2370)]
    (is (= result [2 1]))))

(deftest test-create_blank_tile_tiff
  (with-redefs [gdal/create_geotiff (fn [a b c d e f g h i j k] {:name a :values b :ulx c :uly d :projection e :xdim f :ydim g :xoff h :yoff i})]
    (let [result (raster/create_blank_tile_tiff "foo.tif" 111 222 "wkt_proj" 1 "primary-landcover")]
      (is (= (sort '(:name :values :ulx :uly :projection :xdim :ydim :xoff :yoff)) (sort (keys result))))
      (is (= (count (:values result)) 25000000))
      (is (= 0 (first (:values result)))))))

(deftest test-add_chip_to_tile
  (with-redefs [gdal/update_geotiff (fn [a b c d] {:name a :values b :xoff c :yoff d})]
    (let [result (raster/add_chip_to_tile "foo.tif" [1 2 3] 300 600 330 570)]
      (is (= result {:name "foo.tif" :values [1 2 3] :xoff 1 :yoff 1})))))

(deftest test-nlcd_filter
  (with-redefs [chipmunk/nlcd_filters (fn [a b] {:mask [0 1 0 1] :values [1 2 3 4]})]
    (let [change_result (raster/nlcd_filter [5 6 7 8] "time-since-change" 111 222)
          cover_result  (raster/nlcd_filter [5 6 7 8] "primary-landcover" 111 222)]
      (is (= change_result [0 6 0 8]))
      (is (= cover_result  [1 6 3 8])))))

(deftest test-create_raster
  (with-redefs [util/get-projection           (fn [] "wkt-proj")
                raster/create_blank_tile_tiff (fn [a b c d e f] true)
                storage/ppath                 (fn [a b c d e] (str a b c d e))
                storage/get_url               (fn [a b] (format "http://aws.com/%s/%s" a b) )
                storage/get_json              (fn [i] {"values" [1 2 3 4]})
                raster/nlcd_filter            (fn [a b c d e] [6 7 8 9])
                raster/add_chip_to_tile       (fn [a b c d e f] true)
                storage/put_tiff              (fn [a b] true)
                io/delete-file                (fn [i] true)
                chipmunk/nlcd_filters         (fn [a b] {:mask [0 0] :values [1 2]})]
    (let [input {:date "2007-01-01" :tile "001002" :tilex 111111 :tiley 222222 
                 :chips [{:cx 1 :cy 2 :value 3} {:cx 1 :cy 3 :value 4}] 
                 :product "change"}
          result (raster/create_raster input)]
      (is (string/includes? (first result) "http://aws.com/bar/raster/2007//001/002/change/LCMAP--001002-2007-")) ; (config :storage-destination) defined in project.clj for test
      (is (string/includes? (first result) "-SCTIME.tif"))
      (is (= (count result) 5)))))

(deftest test-create_raster_exception
  (with-redefs [util/get-projection           (fn [] "wkt-proj")
                raster/create_blank_tile_tiff (fn [a b c d e f] true)
                storage/ppath                 (fn [a b c d e] (str a b c d e))
                storage/get_url               (fn [a b] (format "http://aws.com/%s/%s" a b) )
                storage/get_json              (fn [i] {"values" [1 2 3 4]})
                raster/nlcd_filter            (fn [a b c d e] [6 7 8 9])
                raster/add_chip_to_tile       (fn [a b c d e f] true)
                chipmunk/nlcd_filters         (fn [a b] {:mask [0 0] :values [1 2]})
                storage/put_tiff              (fn [a b] (throw (Exception. "exception putting tiff")))
                io/delete-file                (fn [i] true)]
    (let [input {:date "2007-01-01" :tile "001002" :tilex 111111 :tiley 222222 
                 :chips [{:cx 1 :cy 2 :value 3} {:cx 1 :cy 3 :value 4}] 
                 :product "change"}]
      (is (thrown-with-msg? Exception #"exception putting tiff" (raster/create_raster input))))))

(deftest map-details-test
  (with-redefs [config {:region "CU" :ccd_ver "V01"}
                storage/get_url (fn [a b] (str a "/" b))]
    (let [tileid "123456"
          product "change"
          date "2007-07-01"
          product_info (raster/get-products product)
          result (raster/map-details tileid (first product_info) date product)]
      (is (= (keys result) '(:name :prefix :url :data-type :data-product)))
      (is (= (:prefix result) "raster/2007/CU/123/456/change"))
      (is (string/includes? (:name result) "LCMAP-CU-123456-2007-"))
      (is (string/includes? (:name result) "-V01-SCTIME.tif")))))

(deftest product_details_test
  (is (= (keys raster/product_details) 
         '("primary-landcover" "time-of-change" "primary-confidence" "annual-change" "length-of-segment" 
           "curve-fit" "secondary-confidence" "magnitude-of-change" "time-since-change" "secondary-landcover")))
  
  (is (= (map :abbr (vals raster/product_details)) 
         '("LCPRI" "SCTIME" "LCPCONF" "LCACHG" "SCSTAB" "SCMQA" "LCSCONF" "SCMAG" "SCLAST" "LCSEC")))

  (is (= 10 (count raster/product_details))))

