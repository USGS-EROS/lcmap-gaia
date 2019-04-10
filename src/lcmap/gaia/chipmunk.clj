(ns lcmap.gaia.chipmunk
  (:import java.util.Base64)
  (:require [org.httpkit.client    :as http]
            [cheshire.core         :as json]
            [clojure.tools.logging :as log]
            [lcmap.gaia.config     :refer [config]]))

(defn chips_url
  "Return the url for requesting chip data for a given ubid, x, and y"
  ([host ubid x y acquired]
   (str host "/chips?ubid=" ubid "&x=" x "&y=" y "&acquired=" acquired))
  ([ubid x y]
   (chips_url (:chipmunk_host config) ubid x y (:chipmunk_acquired config))))

(defn nlcd
  "Return decoded NLCD chip values for the given X and Y coordinates"
  [x y]
  (let [ubid "AUX_NLCD"
        url (chips_url ubid x y)
        response @(http/get url)
        encoded (-> (:body response) (json/decode) (first) (get "data"))
        decoded (.decode (Base64/getDecoder) encoded)]
    (mapv int decoded)))



