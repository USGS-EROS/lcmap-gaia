(ns lcmap.gaia.server-test
  (:require [clojure.test :refer :all]
            [clojure.string :as string]
            [lcmap.gaia.ccdc   :as ccdc]
            [lcmap.gaia.file   :as file]
            [lcmap.gaia.server :as server]))

(def foo_data (file/read-json "resources/cx-2115585_cy3119805_segment.json"))

(deftest get-product-test-invalid
  (let [response (server/get-product "foo" 111 222 "1990-07-01" {:headers {"accept" "bar"}})]
    (is (= (:status response) 406))
    (is (string/includes? (:body response) "please define a valid Accept header"))))

(deftest get-product-test-valid
  (with-redefs [ccdc/results (fn [x y] foo_data)]
    (let [response (server/get-product "time-since-change" 111 222 "1990-07-01" {:headers {"accept" "application/json"}})
          response_body_keys (set (keys (:body response)))]
      (is (= response_body_keys (set '("x" "y" "values"))))
      (is (= 200 (:status response)))
      (is (= 10000 (count (get (:body response) "values")))))))

(deftest healthy-test
  (let [response (server/healthy {})]
    (is (= (:body response) "OK"))
    (is (= (:status response) 200))))
