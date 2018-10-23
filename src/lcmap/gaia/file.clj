(ns lcmap.gaia.file
  (:gen-class)
  (:require [cheshire.core :as json]))

(defn read-json
  "Returns a lazy sequence"
  [infile]
  (json/parse-stream (clojure.java.io/reader infile)))
