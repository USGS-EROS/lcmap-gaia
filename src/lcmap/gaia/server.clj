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

(defn get-configuration
  [request]
  {:status 200 :body config})

(defn raster-gen
  [{:keys [body] :as req}]
  (log/infof "Received /raster request with params: %s" (dissoc body :chips))
  (try
    (let [map_path (raster/create_geotiff body)]
      {:status 200 :body (assoc (dissoc body :chips) :map_name (:name map_path) :map_prefix (:prefix map_path) :map_url (:url map_path))})
    (catch Exception e
      (log/errorf "Exception in server/raster-gen ! args: %s - message: %s - data: %s - stacktrace: %s" body (.getMessage e) (ex-data e) (-> e stacktrace/print-stack-trace with-out-str))
      {:status 500 :body (assoc body :error (str "problem processing /raster request: " (.getMessage e))
                                     :data (ex-data e))})))

(defn raster-fetch
  [{:keys [body] :as req}]
  true)

(defn product-gen
  [{:keys [body] :as req}]
  (try
    (let [results (products/generate body)
          failures (:failures results)]

      (if (true? (empty? failures))
        {:status 200 :body (dissoc results :failures)}
        {:status 400 :body results}))
    (catch Exception e
      (log/errorf "Exception in server/product-gen ! args: %s - message: %s - data: %s stacktrace:  %s" 
                  body (.getMessage e) (ex-data e) (-> e stacktrace/print-stack-trace with-out-str))
      {:status 500 :body (assoc body :error (str "problem processing /product request: " (.getMessage e))
                                     :data (ex-data e))})))

(defn product-fetch
  [{:keys [body] :as req}]
  true)

(defn healthy
  "Handler for checking application health"
  [request]
  (log/debug "GET health")
  {:status 200 :body {:healthy true}})

(compojure/defroutes routes
  (compojure/context "/" request
    (route/resources "/")
    (compojure/GET   "/healthy" [] (healthy request))
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
