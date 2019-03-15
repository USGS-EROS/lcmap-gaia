(ns lcmap.gaia.server-test
  (:require [clojure.test :refer :all]
            [clojure.string :as string]
            [lcmap.gaia.nemo   :as nemo]
            [lcmap.gaia.server :as server]
            [lcmap.gaia.test-resources :as tr]))

;; (deftest get-product-test-invalid
;;   (let [response (server/get-product "foo" "111" "222" "1990-07-01" {:headers {"accept" "bar"}})]
;;     (is (= (:status response) 406))
;;     (is (string/includes? (:body response) "please define a valid Accept header"))))

;; (deftest get-product-test-valid
;;   (with-redefs [nemo/results (fn [x y] {:segments tr/segments_json :predictions tr/predictions_json})]
;;     (let [response (server/get-product "time-since-change" "111" "222" "1990-07-01" {:headers {"accept" "application/json"}})
;;           response_body_keys (set (keys (:body response)))]
;;       (is (= response_body_keys (set '("x" "y" "values"))))
;;       (is (= 200 (:status response)))
;;       (is (= 10000 (count (get (:body response) "values")))))))

;; (deftest get-products-test
;;   (let [response (server/get-products {:foo "bar"})]
;;     (is (= response 
;;           {:status 200, :body ["annual-change" "curve-fit" "length-of-segment" "magnitude-of-change" "primary-landcover" 
;;                                "primary-landcover-confidence" "secondary-landcover" "secondary-landcover-confidence" 
;;                                "time-of-change" "time-since-change"]}))))

(deftest healthy-test
  (let [response (server/healthy {})]
    (is (= (:body response) {"message" "OK"}))
    (is (= (:status response) 200))))
