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

(defmulti get-product
  (fn [_p _x _y _q request] 
    (log/debugf "GET product request: \nproduct - %s \nx - %s y - %s
                 \nquery day %s \nheaders %s" _p _x _y _q (:headers request))
    (-> request (:headers) (get "accept"))))

(defmethod get-product :default
  [product_type x y query_day request]
  ; https://developer.mozilla.org/en-US/docs/Web/HTTP/Status/406
  {:status 406 :body "please define a valid Accept header of either 'application/json' or..."})

(defmethod get-product "application/json"
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

(defn get-products
  [request]
  {:status 200 :body ["annual-change" "curve-fit" "length-of-segment"  "magnitude-of-change" 
                      "primary-landcover" "primary-landcover-confidence" "secondary-landcover"  
                      "secondary-landcover-confidence" "time-of-change" "time-since-change"]})

(defn get-configuration
  [request]
  {:status 200 :body config})

(defn healthy
  "Hello Gaia"
  [request]
  {:status 200 :body {"message" "OK"}})

(defn product-maps
  [{:keys [body] :as req}]
  (let [{:keys [date tile tilex tiley chips product]} body]
    (log/infof "/maps request -- tile: %s | tilex: %s | tiley: %s | chip count: %s | date: %s | product: %s" 
               tile tilex tiley (count chips) date product)

    {:status 200 :body {"message" product}}
    )
)

(defn persist-product
  [product chipx chipy data query_day]
  (try
    (let [segments (:segments data)
          predictions (:predictions data)
          values (products/data segments predictions product query_day) 
          out_name (util/product-output-name product chipx chipy query_day)
          out_data (json/encode {"x" chipx "y" chipy "values" values})]

      (storage/save_json out_name out_data)
      {:chipx chipx :chipy chipy :date query_day :status "success"})

    (catch Exception e (log/errorf "Exception persist-product: %s" e)
           {:chipx chipx :chipy chipy :date query_day :status "fail"})))

(defn persist-product2
  [product chipxy dates]
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
                  out_data (json/encode {"x" chipx "y" chipy "values" values})]]

      (log/infof "storing chip: %s" chipxy)
      (storage/save_json out_name out_data)
      {:chipx chipx :chipy chipy :date query_day :status "success"})

    (catch Exception e (log/errorf "Exception persist-product: %s" e)
           {:chipxy chipxy :dates dates :status "fail"})))


(defn products
  [{:keys [body] :as req}]
  (log/infof "/products PMAP request body: %s" body)
  (let [{:keys [dates chipx chipy product]} body
        input (nemo/results chipx chipy)
        persist #(persist-product product chipx chipy input %)
        results (pmap persist dates)
        failures (filter (fn [i] (= "fail" (:status i))) results)]

    (if (true? (empty? failures))
      {:status 200 :body {:product product :chipx chipx :chipy chipy}}
      {:status 400 :body {:failed failures}}
      )))

(defn products2
  [{:keys [body] :as req}]
  (log/infof "/products2 request body: %s" body)
  (let [{:keys [dates chips product]} body
        ;input (nemo/results chipx chipy)
        persist #(persist-product2 product % dates)
        results (pmap persist chips)
        failures (filter (fn [i] (= "fail" (:status i))) results)]

    (if (true? (empty? failures))
      {:status 200 :body {:product product :chips chips}}
      {:status 400 :body {:failed failures}}
      )))

(compojure/defroutes routes
  (compojure/context "/" request
                     (route/resources "/")
                     (compojure/GET "/" [] (healthy request))
                     (compojure/GET "/available_products" [] (get-products request))
                     (compojure/GET "/configuration" [] (get-configuration request))
                     (compojure/GET "/product" [product_type x y query_day] (get-product product_type x y query_day request))
                     (compojure/POST "/products" [] (products2 request))
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
