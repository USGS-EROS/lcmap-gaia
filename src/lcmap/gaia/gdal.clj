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
  [pixel_array output_name]
  (let [tif_driver  (gdal/GetDriverByName "GTiff")
        tif_dataset (.Create tif_driver output_name 100 100)
        ;tif_band    (.GetRasterBand tif_datast 1)
        ]
    (.SetGeoTransform tif_dataset ) ;(XULCorner,Cellsize,0,YULCorner,0,-Cellsize)
                                    ; chipx, 30, 0, chipy, 0, -30

    ; SetGeoTransform on dataset
    ; SetProjection on dataset (wkt)
    ; GetRasterBand
    ; .WriteRaster band xoff yoff xsize ysize array

   )



)


;; import org.gdal.gdal Driver
;; import org.gdal.gdal Dataset

;; create a dataset with Driver/Create
;; add a layer with Dataset/CreateLayer

;; import 'org.gdal.gdal gdal
;; (def tiff_driver (gdal/GetDriverByName "GTiff"))
;; (def tiff_dataset (.Create tiff_driver "foo.tif" 100 100))




;; unverified
;; (def band (.GetRasterBand tiff_dataset 1))
;; (.WriteRaster band x x x x x)

