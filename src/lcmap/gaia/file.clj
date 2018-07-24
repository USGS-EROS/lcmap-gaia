(ns lcmap.gaia.file
  (:require [cheshire.core :as json]))

(defn read-json
  "Returns a lazy sequence"
  [infile]
  (json/parse-stream (clojure.java.io/reader infile)))
