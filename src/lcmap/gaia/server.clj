(ns lcmap.gaia.server
  (:require [compojure.core :as compojure]
            [compojure.route :as route]
            [ring.middleware.json :as ring-json]
            [ring.middleware.keyword-params :as ring-keyword-params]
            [ring.middleware.defaults :as ring-defaults]
            [ring.util.response :as ring-response]
            [org.httpkit.client :as http]
            [org.httpkit.server :as http-server]
            [lcmap.gaia.file :as file]
            [lcmap.gaia.products :as products]
            [lcmap.gaia.util :as util])
  (:import (java.io FileInputStream)))

(def foo_data (file/read-json "resources/y3161805_x-2115585_nodates.json"))
(def query_day "2006-07-01")

(defmulti get-product
  (fn [_ request] (-> request (:headers) (get "accept"))))

(defmethod get-product :default
  [product_type request]
  {:status 200 :body [(str "please define a valid Accept header of either 'application/json' or...") ]})

(defmethod get-product "application/json"
  [product_type request]
  (let [input foo_data
        product_fn (-> (str "lcmap.gaia.products/" product_type) (symbol) (resolve))
        product_values (products/data input product_fn query_day)
        chipx (get (first input) "chipx")
        chipy (get (first input) "chipy")]
    {:status 200 :body {"chipx" chipx "chipy" chipy "values" product_values}}))

(defn healthy
  "Hello Gaia"
  [request]
  {:status 200 :body ["OK"]})

(compojure/defroutes routes
  (compojure/context "/" request
                     (route/resources "/")
                     (compojure/GET "/" [] (healthy request))
                     (compojure/GET "/product/:product_type" [product_type] (get-product product_type request))))

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
 (http-server/run-server app {:port 9876}))
