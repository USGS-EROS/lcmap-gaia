(ns lcmap.gaia.config-test
  (:require [clojure.test :refer :all]
            [lcmap.gaia.config :as config]
            [environ.core :as environ]))


;; (:storage-destination config) defined in project.clj for test
(deftest storage_destination_defined_test
    (let [cfg config/config]
      (is (= "bar" (:storage-destination cfg)))
      (is (= "foo" (:storage-bucket cfg)))))
