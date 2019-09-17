(ns lcmap.gaia.nemo-test
  (:use org.httpkit.fake)
  (:require [clojure.test      :refer :all]
            [lcmap.gaia.nemo   :as nemo]
            [lcmap.gaia.config :as config]
            [org.httpkit.client :as http]))

(def faux_config {:nemo_host "http://nemo-host.org" :segments_path "/segments" :predictions_path "/predictions"})
(def faux_nemo_segments_url "http://nemo-host.org/segments?cx=999&cy=666")
(def faux_nemo_predictions_url "http://nemo-host.org/predictions?cx=999&cy=666")

(deftest test-results_url
  (with-redefs [config/config faux_config]
    (let [url (nemo/results-url 999 666 (:segments_path faux_config))]
      (is (= url faux_nemo_segments_url)))))

(deftest test-results-good
  (with-fake-http [{:url faux_nemo_segments_url :method :get} {:status 200 :body "{\"xxx\":\"yyy\"}"}
                   {:url faux_nemo_predictions_url :method :get} {:status 200 :body "{\"xxx\":\"yyy\"}"}]
    (with-redefs [config/config faux_config]
      (let [segments (nemo/segments 999 666)
            predictions (nemo/predictions 999 666)]
        (is (= segments {"xxx" "yyy"}))))))

(deftest test-results-bad
  (with-fake-http [{:url faux_nemo_segments_url :method :get} {:status 500 :body nil}
                   {:url faux_nemo_predictions_url :method :get} {:status 404 :body "err"}]
    (with-redefs [config/config faux_config]
      (is (thrown-with-msg? Exception #"Error requesting segments data from Nemo" (nemo/segments 999 666))))))
