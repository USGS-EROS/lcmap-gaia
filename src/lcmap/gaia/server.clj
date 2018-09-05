(ns lcmap.gaia.server
  (:require [compojure.core :as compojure]
            [compojure.route :as route]
            [ring.middleware.json :as ring-json]
            [ring.middleware.keyword-params :as ring-keyword-params]
            [ring.middleware.defaults :as ring-defaults]
            [org.httpkit.client :as http]
            [org.httpkit.server :as http-server]))

(defn healthy
  "Hello Gaia"
  [request]
  {:status 200 :body ["OK"]})

(compojure/defroutes routes
  (compojure/context "/" request
                     (route/resources "/")
                     (compojure/GET "/" [] (healthy request))))

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
