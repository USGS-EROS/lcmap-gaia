(ns lcmap.gaia.chipmunk
  (:require [org.httpkit.client    :as http]
            [cheshire.core         :as json]
            [clojure.tools.logging :as log]
            [lcmap.gaia.config     :refer [config]]))

(defn results_url
  ([x y host resource]
   (str host "?")
   ; $ http GET "http://lcmap-test.cr.usgs.gov/aux_cu_v01/chips?ubid=AUX_NLCD&x=2081415&y=2423805&acquired=1999-01-01/2002-01-01"
   )
  ([x y resource]
   (results_url x y (:chipmunk_host config) resource)))
