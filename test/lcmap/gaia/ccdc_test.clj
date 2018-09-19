(ns lcmap.gaia.ccdc-test
  (:require [clojure.test      :refer :all]
            [lcmap.gaia.ccdc   :as ccdc]
            [lcmap.gaia.config :as config]
            [org.httpkit.fake  :refer [with-fake-http]]))

(def faux_config {:nemo_host "http://nemo-host.org" :nemo_resource "/foo"})
(def faux_nemo_url "http://nemo-host.org/foo?chipx=999&chipy=666")

(deftest test-results_url
  (with-redefs [config/config faux_config]
    (let [url (ccdc/results_url 999 666)]
      (is (= url faux_nemo_url)))))

(deftest test-results
  (with-fake-http [{:url faux_nemo_url :method :get} {:status 200 :body "[\"xxx\"]"}]
    (with-redefs [config/config faux_config]
      (let [results (ccdc/results 999 666)]
        (is (= 200 (:status results)))
        (is (= "xxx" (:body results)))))))
