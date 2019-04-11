(ns lcmap.gaia.server-test
  (:require [clojure.test :refer :all]
            [clojure.string :as string]
            [lcmap.gaia.nemo   :as nemo]
            [lcmap.gaia.server :as server]
            [lcmap.gaia.test-resources :as tr]))

(deftest healthy-test
  (let [response (server/healthy {})]
    (is (= (:body response) {"message" "OK"}))
    (is (= (:status response) 200))))
