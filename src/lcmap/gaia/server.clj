(ns lcmap.gaia.server
  (:require [compojure.core :as compojure]
            [compojure.route :as route]
            [ring.middleware.json :as ring-json]
            [ring.middleware.keyword-params :as ring-keyword-params]
            [ring.middleware.defaults :as ring-defaults]
            [ring.util.response :as ring-response]
            [org.httpkit.client :as http]
            [org.httpkit.server :as http-server]
            [lcmap.gaia.ccdc :as ccdc]
            [lcmap.gaia.file :as file]
            [lcmap.gaia.products :as products]
            [lcmap.gaia.util :as util]))


(defmulti get-product
  (fn [_p _q _x _y request] (-> request (:headers) (get "accept"))))

(defmethod get-product :default
  [product_type x y query_day request]
  ; https://developer.mozilla.org/en-US/docs/Web/HTTP/Status/406
  {:status 406 :body "please define a valid Accept header of either 'application/json' or..."})

(defmethod get-product "application/json"
  [product_type x y query_day request]
  (let [input (ccdc/results x y) 
        product_values (products/data input product_type query_day)]
    {:status 200 :body {"x" x "y" y "values" product_values}}))

(defn get-products
  [request]
  {:status 200 :body ["curve-fit" "time-of-change" "time-since-change" 
                      "magnitude-of-change" "length-of-segment"]})

(defn healthy
  "Hello Gaia"
  [request]
  {:status 200 :body "OK"})

(compojure/defroutes routes
  (compojure/context "/" request
                     (route/resources "/")
                     (compojure/GET "/" [] (healthy request))
                     (compojure/GET "/available-products" [] (get-products))
                     (compojure/GET "/product/:product_type/:x/:y/:query_day" [product_type x y query_day] (get-product product_type x y query_day request))))

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
