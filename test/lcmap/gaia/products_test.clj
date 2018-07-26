(ns lcmap.gaia.products-test
  (:require [clojure.test :refer :all]
            [lcmap.gaia.products :as products]
            [lcmap.gaia.util     :as util]))

(deftest time-of-change-positive-test
  (let [queryday "1990-04-01"
        breakord (-> "1990-03-01" (util/to-javatime) (util/javatime-to-ordinal))
        model {"bday" breakord "chprob" 1.0}]
    (is (= 111 (products/time-of-change model queryday)))
    )
)

;; (deftest time-of-change-negative-test)
