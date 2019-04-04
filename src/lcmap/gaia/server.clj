(ns lcmap.gaia.server
  (:require [clojure.tools.logging :as log]
            [clojure.stacktrace :as stacktrace]
            [compojure.core :as compojure]
            [compojure.route :as route]
            [ring.middleware.json :as ring-json]
            [ring.middleware.keyword-params :as ring-keyword-params]
            [ring.middleware.defaults :as ring-defaults]
            [ring.util.response :as ring-response]
            [org.httpkit.client :as http]
            [org.httpkit.server :as http-server]
            [lcmap.gaia.nemo :as nemo]
            [lcmap.gaia.file :as file]
            [lcmap.gaia.products :as products]
            [lcmap.gaia.util :as util]
            [lcmap.gaia.config :refer [config]]
            [lcmap.gaia.raster :as raster]
            [lcmap.gaia.storage :as storage]
            [cheshire.core :as json]))

(defn healthy
  "Hello Gaia"
  [request]
  {:status 200 :body {"message" "OK"}})

(defn available-products
  [request]
  {:status 200 :body ["annual-change" "curve-fit" "length-of-segment"  "magnitude-of-change" 
                      "primary-landcover" "primary-landcover-confidence" "secondary-landcover"  
                      "secondary-landcover-confidence" "time-of-change" "time-since-change"]})

(defn get-configuration
  [request]
  {:status 200 :body config})

(defn raster-gen
  [{:keys [body] :as req}]
  (try
    (let [map_path (raster/create-geotiff body)]
      {:status 200 :body (assoc (dissoc body :chips) :map_name (:name map_path))})
    (catch Exception e
      (log/errorf "Exception in product-maps: %s" (ex-data e))
      {:status 500 :body {:error (str "problem processing /raster request: " (ex-data e)) 
                          :body-minus-chips (:dissoc body :chips)}})))

(defn raster-fetch
  [{:keys [body] :as req}]
  true)

(defn product-gen
  [{{dates :dates cx :cx cy :cy product :product tile :tile} :body :as all}]
  (try
    (let [results (products/generation all)
          failures (:failures results)]

      (if (true? (empty? failures))
        {:status 200 :body (dissoc results :failures)}
        {:status 400 :body results}))
    (catch Exception e
      (log/errorf "Exception in product-gen: %s" (-> e stacktrace/print-stack-trace with-out-str))
      {:status 500 :body (assoc (:body all) :error (str "problem processing /product request: " (.getMessage e)))})))

(defn product-fetch
  [{:keys [body] :as req}]
  true)

(compojure/defroutes routes
  (compojure/context "/" request
    (route/resources "/")
    (compojure/GET   "/" [] (healthy request))
    (compojure/GET   "/available_products" [] (available-products request))
    (compojure/GET   "/configuration" [] (get-configuration request))
    (compojure/POST  "/product"  [] (product-gen   request))
    (compojure/GET   "/product"  [] (product-fetch request))
    (compojure/POST  "/raster"   [] (raster-gen    request))
    (compojure/GET   "/raster"   [] (raster-fetch  request))))

(defn response-handler
  [routes]
  (-> routes
      (ring-json/wrap-json-body {:keywords? true})
      (ring-json/wrap-json-response)
      (ring-defaults/wrap-defaults ring-defaults/api-defaults)
      (ring-keyword-params/wrap-keyword-params)))

(def app (response-handler routes))

(defn run-server
  []
 (http-server/run-server app {:port (:http_port config)}))
