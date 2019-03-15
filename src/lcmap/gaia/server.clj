(ns lcmap.gaia.server
  (:require [clojure.tools.logging :as log]
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

(defn product-maps
  [{:keys [body] :as req}]
  (try
    (let [{:keys [date tile tilex tiley chips product]} body
          map_name (products/map-name tile product date)
          projection (util/get-projection)]
      (raster/create_blank_tile_tiff map_name tilex tiley projection)
      (doseq [chip chips
              :let [cx (:cx chip)
                    cy (:cy chip)
                    chip_name (util/product-output-name product cx cy date)
                    chip_data (storage/get_json tile chip_name)]]
        (log/debugf "adding %s to tile: %s" chip_name tile)
        (if chip_data
          (raster/add_chip_to_tile map_name (get chip_data "values") tilex tiley cx cy)
          (log/debugf "no data to add to tile %s at cx: %s | cy: %s" tile cx cy)))
      (storage/put_tiff tile map_name map_name)
      (log/infof "done adding data to map_name: %s" map_name)
      {:status 200 :body (assoc (dissoc body :chips) :map_name map_name)})
    (catch Exception e
      (log/errorf "Exception in product-maps: %s" (ex-data e))
      {:status 500 :body {:error (str "problem processing /maps request: " (ex-data e)) :body-minus-chips (:dissoc body :chips)}})))

(defn persist-product
  [product cx cy query_day tile data]
  (try
    (log/infof "working on cx: %s  cy: %s  date: " cx cy query_day)
    (let [values (products/data (:segments data) (:predictions data) product query_day) 
          out_name (util/product-output-name product cx cy query_day)
          out_data {"x" cx "y" cy "values" values}]

      (log/infof "storing cx: %s  cy: %s  date: %s" cx cy query_day)
      (storage/put_json tile out_name out_data)
      {:chipx cx :chipy cy :date query_day :status "success"})

    (catch Exception e (log/errorf "Exception persist-product: %s" e)
           {:cx cx :cy cy :date query_day :product product :status "fail"})))

(defn products
  [{:keys [body] :as req}]
  (try
    (let [{:keys [dates cx cy product tile]} body
          __ (storage/create_bucket tile)
          data (nemo/results cx cy)
          persist #(persist-product product cx cy % tile data)
          results (pmap persist dates)
          failures (filter (fn [i] (= "fail" (:status i))) results)]

      (if (true? (empty? failures))
        {:status 200 :body {:product product :cx cx :cy cy :dates dates}}
        {:status 400 :body {:failed failures}}))
    (catch Exception e
      (log/errorf "Exception in products: %s" (ex-data e))
      {:status 500 :body {:error (str "problem processing /products request: " (ex-data e)) :request-body body}})))

(compojure/defroutes routes
  (compojure/context "/" request
    (route/resources "/")
    (compojure/GET "/" [] (healthy request))
    (compojure/GET "/available_products" [] (available-products request))
    (compojure/GET "/configuration" [] (get-configuration request))
    (compojure/POST "/products" [] (products request))
    (compojure/POST "/maps" [] (product-maps request))))

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
