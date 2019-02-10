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
            [cheshire.core :as json]
            ))

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

(defn get-product
  [product_type x y query_day request]
  (try
    (let [input (nemo/results x y) 
          product_values (products/data (:segments input) (:predictions input) product_type query_day)]
      {:status 200 :body {"x" (read-string x) "y" (read-string y) "values" product_values}})
    (catch Exception e
      (cond 
        (= :data-failure (-> e ex-data :cause))
          (do (log/errorf "data request problem: %s" (ex-data e))
              {:status 500 :body "Unable to retrieve input data"})
        (= :validation-failure (-> e ex-data :cause))
          (do (log/errorf "input validation problem: %s" (ex-data e))
              {:status 400 :body "Invalid input data"})
        :else
          (do (log/errorf "Exception encountered handling get-product request: %s" (ex-data e))
              {:status 500 :body "Error handling request"})))))

(defn product-maps
  [{:keys [body] :as req}]
  (let [{:keys [date tile tilex tiley chips product]} body
        tile_name (str tile "_" product "_" date ".tif")
        projection (util/get-projection)]
    (log/infof "/maps request -- tile: %s | tilex: %s | tiley: %s | chip count: %s | date: %s | product: %s" 
               tile tilex tiley (count chips) date product)
    ; create blank tile tiff. in local tmp dir?
    (raster/create_blank_tile_tiff tile_name tilex tiley projection)
    ; loop through coordinate pairs, calculating chip tiff name, checking for existence, and grabbing values
    ; when available, update tiff in local tmp dir
    (doseq [chip chips
            :let [cx (:cx chip)
                  cy (:cy chip)
                  chip_name (util/product-output-name product cx cy date)
                  chip_data (storage/get_json tile name)
                  ]]
      (raster/add_chip_to_tile tile_name chip_data tilex tiley cx cy)
      )

    ; when complete, store tiff in object store
    ; remove tiff from local tmp dir



    {:status 200 :body {"message" product}}
    )
)

(defn persist-product
  [product chipxy dates tile]
  (try
    (log/infof "working on chip: %s" chipxy)
    (doseq [query_day dates
            :let [chipx (:cx chipxy)
                  chipy (:cy chipxy)
                  data (nemo/results chipx chipy)
                  segments (:segments data)
                  predictions (:predictions data)
                  values (products/data segments predictions product query_day) 
                  out_name (util/product-output-name product chipx chipy query_day)
                  out_data {"x" chipx "y" chipy "values" values}]]

      (log/infof "storing chip: %s" chipxy)
      (storage/put_json tile out_name out_data)
      {:chipx chipx :chipy chipy :date query_day :status "success"})

    (catch Exception e (log/errorf "Exception persist-product: %s" e)
           {:chipxy chipxy :dates dates :status "fail"})))

(defn products
  [{:keys [body] :as req}]
  (log/infof "/products2 request body: %s" body)
  (let [{:keys [dates chips product tile]} body
        __ (storage/create_bucket tile)
        persist #(persist-product product % dates tile)
        results (pmap persist chips)
        failures (filter (fn [i] (= "fail" (:status i))) results)]

    (if (true? (empty? failures))
      {:status 200 :body {:product product :chips chips}}
      {:status 400 :body {:failed failures}})))

(compojure/defroutes routes
  (compojure/context "/" request
    (route/resources "/")
    (compojure/GET "/" [] (healthy request))
    (compojure/GET "/available_products" [] (available-products request))
    (compojure/GET "/configuration" [] (get-configuration request))
    (compojure/GET "/product" [product_type x y query_day] (get-product product_type x y query_day request))
    (compojure/POST "/products" [] (products request))
    (compojure/POST "/maps" [] (product-maps request))))

(defn response-handler
  [routes]
  (-> routes
      (ring-json/wrap-json-body {:keywords? true})
      (ring-json/wrap-json-response)
      (ring-defaults/wrap-defaults ring-defaults/api-defaults)
      (ring-keyword-params/wrap-keyword-params))
)
(def app (response-handler routes))

(defn run-server
  []
 (http-server/run-server app {:port (:http_port config)}))
