(ns lcmap.gaia.nemo
  (:require [org.httpkit.client    :as http]
            [cheshire.core         :as json]
            [clojure.tools.logging :as log]
            [lcmap.gaia.config     :refer [config]]))

(defn results_url
  ([x y host path]
   (str host path "?cx=" x "&cy=" y))
  ([x y path]
   (results_url x y (:nemo_host config) path)))

(defn parse_body
  [http_response]
  (json/parse-string (:body http_response)))

(defn segments
  [x y]
  (let [url (results_url x y (:segments_path config))
        response @(http/get url)]
    (if (= 200 (:status response))
      (parse_body response)
      (throw (ex-info "Nemo Error - non-200 /segment response" 
                      {:type :nemo-request-failure :cause :data-failure 
                       :status (:status response) :url url})))))

(defn predictions
  [x y]
  (let [url (results_url x y (:predictions_path config))
        response @(http/get url)]
    (if (= 200 (:status response))
      (parse_body response)
      (throw (ex-info "Nemo Error - non-200 /annual_prediction response" 
                      {:type :nemo-request-failure :cause :data-failure 
                       :status (:status response) :url url})))))

