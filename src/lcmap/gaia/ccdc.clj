(ns lcmap.gaia.ccdc
  (:require [org.httpkit.client    :as http]
            [cheshire.core         :as json]
            [clojure.tools.logging :as log]
            [lcmap.gaia.config     :refer [config]]))

(defn results_url
  [x y]
  (str (:nemo_host config) (:nemo_resource config) "?chipx=" x "&chipy=" y))

(defn results
  [x y]
  (let [response @(http/get (results_url x y))
        status (:status response)]
    (if (= 200 status)
      (json/parse-string (:body response))
      (do (log/errorf "problem retrieving ccdc results. status code: %s \n response %s " status response)
          false))))
