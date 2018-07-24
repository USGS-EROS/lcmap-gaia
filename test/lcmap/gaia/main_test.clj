(ns lcmap.gaia.main-test
  (:require [clojure.test :refer :all]
            [lcmap.gaia.main :as gm]))


(deftest is-quad?
  (is (= 4 (gm/foo 2))))
