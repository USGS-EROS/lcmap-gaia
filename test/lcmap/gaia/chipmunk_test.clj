(ns lcmap.gaia.chipmunk-test
  (:use org.httpkit.fake)
  (:import java.util.Base64)
  (:require [clojure.test        :refer :all]
            [cheshire.core       :as json]
            [lcmap.gaia.chipmunk :as chipmunk]))

(def nlcd-cx 111111)
(def nlcd-cy 222222)
(def nlcd-url (str "http://localhost:5656/chips?ubid=AUX_NLCD&x=" nlcd-cx "&y=" nlcd-cy "&acquired=1999-01-01/2002-01-01"))
(def nlcd-vals [41 52 63])

(defn nlcd_resp []
  (let [byte-vals (byte-array nlcd-vals)
        enc-vals (.encodeToString (Base64/getEncoder) byte-vals)]
    (json/encode [{"data" enc-vals}])))

(deftest test-chips_url
  (let [url (chipmunk/chips_url "AUX_NLCD" nlcd-cx nlcd-cy)]
    (is (= url nlcd-url))))

(deftest test-nlcd
  (with-fake-http [{:url nlcd-url :method :get} {:status 200 :body (nlcd_resp)}]
    (let [response (chipmunk/nlcd nlcd-cx nlcd-cy)]
      (is (= response nlcd-vals)))))


