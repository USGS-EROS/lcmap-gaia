(ns faux.main
  (:gen-class)
  (:require [compojure.core :as compojure]
            [compojure.route :as route]
            [cheshire.core :as json]
            [ring.middleware.json :as ring-json]
            [ring.middleware.keyword-params :as ring-keyword-params]
            [ring.middleware.defaults :as ring-defaults]
            [ring.util.response :as ring-response]
            [org.httpkit.client :as http]
            [org.httpkit.server :as http-server]))

(def predictions_json "cx-2115585_cy3119805_prediction_with_fake_date.json")
(def segments_json "cx-2115585_cy3119805_segment.json")

(defn read-json
  "Returns a lazy sequence"
  [infile]
  (json/parse-stream (clojure.java.io/reader infile)))

(defn send_json
  [dtype]
  (if (= dtype "predictions")
    (read-json predictions_json)
    (read-json segments_json)))

(defn healthy
  [request]
  {:status 200 :body {"message" "OK"}})

(compojure/defroutes routes
  (compojure/context "/" request
                     (route/resources "/")
                     (compojure/GET "/" [] (healthy request))
                     (compojure/GET "/predictions" [] (send_json "predictions"))
                     (compojure/GET "/segments" [] (send_json "segments"))))

(defn response-handler
  [routes]
  (-> routes
      (ring-json/wrap-json-body {:keywords? true})
      (ring-json/wrap-json-response)
      (ring-defaults/wrap-defaults ring-defaults/api-defaults)
      (ring-keyword-params/wrap-keyword-params)))

(def app (response-handler routes))

(defn -main
  []
 (http-server/run-server app {:port 8899}))

