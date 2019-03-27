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

(defn product-maps
  [{:keys [body] :as req}]
  (try
    (let [{:keys [date tile tilex tiley chips product]} body
          map_path (products/map-path tile product date)
          projection (util/get-projection)]
      (raster/create_blank_tile_tiff (:name map_path) tilex tiley projection)
      (log/infof "Received /maps request to produce: %s" (:name map_path))
      (doseq [chip chips
              :let [cx (:cx chip)
                    cy (:cy chip)
                    chip_path (products/product-path product cx cy date)
                    chip_data (storage/get_json chip_path)]]
        (log/debugf "adding %s to tile: %s" (:name chip_path) tile)
        (if chip_data
          (raster/add_chip_to_tile (:name map_path) (get chip_data "values") tilex tiley cx cy)
          (log/debugf "no data to add to tile %s at cx: %s | cy: %s" tile cx cy)))
      (storage/put_tiff map_path (:name map_path))
      (log/infof "done adding data to map_name: %s" (:name map_path))
      {:status 200 :body (assoc (dissoc body :chips) :map_name (:name map_path))})
    (catch Exception e
      (log/errorf "Exception in product-maps: %s" (ex-data e))
      {:status 500 :body {:error (str "problem processing /maps request: " (ex-data e)) :body-minus-chips (:dissoc body :chips)}})))

(defn persist-product
  [product cx cy tile query_day data]
  (try
    (let [values (products/data (:segments data) (:predictions data) product query_day) 
          out_path (products/product-path product cx cy tile query_day)
          out_data {"x" cx "y" cy "values" values}]

      (log/infof "storing : %s" (:name out_path))
      (storage/put_json out_path out_data)
      {:cx cx :cy cy :date query_day :status "success"})

    (catch Exception e (log/errorf "Exception in persist-product - cx: %s  cy: %s  product: %s  date: %s exception-message: %s exception-data: %s" 
                                   cx cy product query_day (.getMessage e) (ex-data e))
           {:cx cx :cy cy :date query_day :product product :status "fail" :message (.getMessage e)})))

(defn products
  [{{dates :dates cx :cx cy :cy product :product tile :tile} :body :as all}]
  (try
    (let [data (nemo/results cx cy)
          persist #(persist-product product cx cy tile % data)
          results (pmap persist dates)
          failures (->> results (filter (fn [i] (= "fail" (:status i)))) (map (fn [i] {(:date i) (:message i)})))
          response_map {:product product :cx cx :cy cy :dates dates}]

      (if (true? (empty? failures))
        {:status 200 :body response_map}
        {:status 400 :body (assoc response_map :failed_dates failures)}))
    (catch Exception e
      (log/errorf "Exception in products: %s" (-> e stacktrace/print-stack-trace with-out-str))
      {:status 500 :body (assoc (:body all) :error (str "problem processing /products request: " (.getMessage e)))})))

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
