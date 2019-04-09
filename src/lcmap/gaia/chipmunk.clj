(ns lcmap.gaia.chipmunk
  (:require [org.httpkit.client    :as http]
            [cheshire.core         :as json]
            [clojure.tools.logging :as log]
            [lcmap.gaia.config     :refer [config]]))

(defn results_url
  ([host resource ubid x y acquired]
   (str host "/" resource "?ubid=" ubid "&x=" x "&y=" y "&acquired=" acquired))
  ([resource ubid x y]
   (results_url (:chipmunk_host config) resource ubid x y (:chipmunk_acquired config))))

(defn nlcd
  [x y]
  (let [url (results_url "chips" "AUX_NLCD" x y)
        response @(http/get url)
        body (first (json/decode (:body response)))]
    (get body "data")))
