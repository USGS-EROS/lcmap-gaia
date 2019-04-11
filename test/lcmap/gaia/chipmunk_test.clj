(ns lcmap.gaia.chipmunk-test
  (:use org.httpkit.fake)
  (:import java.util.Base64)
  (:require [clojure.test        :refer :all]
            [cheshire.core       :as json]
            [lcmap.gaia.chipmunk :as chipmunk]))

(def nlcd-cx 111111)
(def nlcd-cy 222222)
(def nlcd-url (str "http://localhost:5656/chips?ubid=AUX_NLCD&x=" nlcd-cx "&y=" nlcd-cy "&acquired=1999-01-01/2002-01-01"))
(def nlcd-vals [41 52 0])
(def nlcd-mask-vals [1 1 0])

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

(deftest test-binary
  (is (= 1 (chipmunk/binary 55)))
  (is (= 0 (chipmunk/binary 0))))

(deftest test-nlcd_mask
  (with-fake-http [{:url nlcd-url :method :get} {:status 200 :body (nlcd_resp)}]
    (let [response (chipmunk/nlcd_mask nlcd-cx nlcd-cy)]
      (is (= response nlcd-mask-vals))))
  (is (= nlcd-mask-vals (chipmunk/nlcd_mask nlcd-vals))))

(deftest test-nlcd_filters
  (with-fake-http [{:url nlcd-url :method :get} {:status 200 :body (nlcd_resp)}]
    (let [response (chipmunk/nlcd_filters nlcd-cx nlcd-cy)]
      (is (= (:values response) nlcd-vals))
      (is (= (:mask response) nlcd-mask-vals))))



)
