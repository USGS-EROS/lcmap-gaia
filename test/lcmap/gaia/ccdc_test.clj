(ns lcmap.gaia.ccdc-test
  (:use org.httpkit.fake)
  (:require [clojure.test      :refer :all]
            [lcmap.gaia.ccdc   :as ccdc]
            [lcmap.gaia.config :as config]
            [org.httpkit.client :as http]))

(def faux_config {:nemo_host "http://nemo-host.org" :nemo_resource "/foo"})
(def faux_nemo_url "http://nemo-host.org/foo?chipx=999&chipy=666")

(deftest test-results_url
  (with-redefs [config/config faux_config]
    (let [url (ccdc/results_url 999 666)]
      (is (= url faux_nemo_url)))))

(deftest test-results-good
  (with-fake-http [{:url faux_nemo_url :method :get} {:status 200 :body "{\"xxx\":\"yyy\"}"}]
    (with-redefs [config/config faux_config]
      (let [results (ccdc/results 999 666)]
        (println (str "results: " results))
        (is (= {"xxx" "yyy"} results))))))

(deftest test-results-bad
  (with-fake-http [{:url faux_nemo_url :method :get} {:status 404 :body "resource not found"}]
    (with-redefs [config/config faux_config]
      (let [results (ccdc/results 999 666)]
        (is (= false results))))))
