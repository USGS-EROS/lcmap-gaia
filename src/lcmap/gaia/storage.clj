(ns lcmap.gaia.storage
  (:require [lcmap.gaia.config :refer [config]]))


(defn list_buckets
  []
  true)

(defn create_bucket
  [name]
  true)

(defn add_file
  [filename bucket]
  true)

(defn drop_file
  [filename bucket]
  true)

(defn save_json
  [filename data]
  (let [output (str (:output_path config) filename)]
    (spit output data)))
