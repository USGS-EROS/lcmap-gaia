(ns lcmap.gaia.server-test
  (:require [clojure.test :refer :all]
            [clojure.string :as string]
            [lcmap.gaia.server :as server]))

(deftest get-product-test-invalid
  (let [response (server/get-product "foo" 111 222 "1990-07-01" {:headers {"accept" "bar"}})]
    (is (= (:status response) 406))
    (is (string/includes? (:body response) "please define a valid Accept header"))))

(deftest get-product-test-json
  (let [response (server/get-product "time-since-change" 111 222 "1990-07-01" {:headers {"accept" "application/json"}})
        response_body_keys (set (keys (:body response)))]
    (is (= response_body_keys (set '("x" "y" "values"))))))

(deftest healthy-test
  (let [response (server/healthy {})]
    (is (= (:body response) "OK"))
    (is (= (:status response) 200))))
