(ns lcmap.gaia.gdal
  (:require [mount.core :as mount]
            [lcmap.gaia.util :as util])
  (:import [org.gdal.gdal gdal]
           [org.gdal.gdal Driver]
           [org.gdal.gdal Dataset]
           [org.gdal.gdalconst gdalconst]))

;; init and state constructs blatantly ripped off from the
;; USGS-EROS/lcmap-chipmunk project on GitHub, created by
;; Jon Morton https://github.com/jmorton
;; 
;; ## Init
;;
;; This makes it easier to use Java GDAL libraries without
;; having to set environment variables. These are typical
;; install locations of GDAL libs on CentOS and Ubuntu.
;;
;; Before GDAL can open files, drivers must be registered.
;; Selective registration is more tedious and error prone,
;; so we just register all drivers.
;;
;; If anything goes wrong, a helpful string is printed to
;; stdout (not a log file).
;;

(defn init
  "Initialize GDAL drivers."
  []
  (try
    (util/amend-usr-path ["/usr/lib/java/gdal" "/usr/lib/jni"])
    (gdal/AllRegister)
    (catch RuntimeException e
      (binding [*out* *err*]
        (println (str "Could not update paths to native libraries. "
                      "You may need to set LD_LIBRARY_PATH to the "
                      "directory containing libgdaljni.so"))))
    (finally
      (import org.gdal.gdal.gdal))))


;; ## State
;;
;; A mount state is defined so that GDAL is initialized like
;; everything else (DB connections, HTTP listeners, etc...)
;;

(mount/defstate gdal-init
  :start (init))

(defn geotiff_from_pixel_array
  [pixel_array output_name ulx uly proj_wkt]
  (let [tif_driver  (gdal/GetDriverByName "GTiff")
        tif_dataset (.Create tif_driver output_name 100 100)
        tif_band    (.GetRasterBand tif_dataset 1)
        transform_array (double-array [ulx 30 0 uly 0 30])] ;(XULCorner,Cellsize,0,YULCorner,0,-Cellsize)
    (.SetGeoTransform tif_dataset transform_array) 
    (.SetProjection tif_dataset proj_wkt)
    (.WriteRaster tif_band 0 0 100 100 (float-array pixel_array)) 
    ; write to disk
    (.delete tif_band)
    (.delete tif_dataset))
  output_name)

