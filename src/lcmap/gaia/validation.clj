(ns lcmap.gaia.validation
  (:require [clojure.spec.alpha :as spec]))

(defn check!
  "Return args if they conform to spec"
  [spec params]
  (or (some->> (spec/explain-data spec params)
               (ex-info "validation error")
               (throw))
      (spec/conform spec params)))
