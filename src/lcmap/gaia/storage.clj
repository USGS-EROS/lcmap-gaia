(ns lcmap.gaia.storage
  (:require [mount.core :as mount]
            [clojure.tools.logging :as log]
            [lcmap.gaia.config :refer [config]]
            [cheshire.core :as json]
            [amazonica.aws.s3 :as s3]
            [java-time :as jt]))

(def bucketname (:storage-bucket config))

(def client-config
  {:access-key (:storage-access-key config)
   :secret-key (:storage-secret-key config)
   :endpoint   (:storage-endpoint config)
   :client-config {:path-style-access-enabled true}})

(defn list_buckets
  ([cfg]
   (let [objects (s3/list-buckets cfg)]
     (map :name objects)))
  ([]
   (list_buckets client-config)))

(defn list_bucket_contents
  ([bucket prefix]
   (let [objects (s3/list-objects-v2 client-config :bucket-name bucket :prefix prefix :delimiter "/")
         summaries (:object-summaries objects)]
     (map :key summaries)))
  ([bucket]
   (let [objects (s3/list-objects-v2 client-config :bucket-name bucket)
         summaries (:object-summaries objects)]
     (map :key summaries)))
  ([]
   (list_bucket_contents bucketname)))

(defn list_bucket_prefixes
  ([bucket prefix]
   (:common-prefixes (s3/list-objects-v2 client-config :bucket-name bucket :prefix prefix :delimiter "/")))
  ([bucket]
   (list_bucket_prefixes bucket ""))
  ([]
   (list_bucket_prefixes bucketname)))

(defn create_bucket
  ([bucket]
   (s3/create-bucket client-config bucket))
  ([]
   (create_bucket bucketname)))

(defn put_json
  ([bucket output_path data]
   (let [encoded_data (json/encode data)
         byte_data (.getBytes encoded_data)
         byte_stream (java.io.ByteArrayInputStream. byte_data)
         metadata {:content-length (count byte_data)}
         keyname (str (:prefix output_path) "/" (:name output_path))]
     (try
       (s3/put-object client-config :bucket-name bucket :key keyname :input-stream byte_stream :metadata metadata)
       true
       (catch Exception e (log/errorf "Error putting data to object store: %s" e)
              false))))
  ([output_path data]
   (put_json bucketname output_path data)))

(defn put_tiff
  ([bucket filepath filelocation]
   (let [javafile (java.io.File. filelocation)
         keyname (str (:prefix filepath) "/" (:name filepath))]
     (try
       (s3/put-object client-config :bucket-name bucket :key keyname  :file javafile)
       true
       (catch Exception e (log/errorf "Error putting data to object store: %s" e)
              false))))
  ([filepath filelocation]
   (put_tiff bucketname filepath filelocation)))

(defn drop_object
  ([bucket filename]
   (s3/delete-object client-config :bucket-name bucket :key filename))
  ([filename]
   (drop_object bucketname filename)))

(defn drop_bucket
  [bucket]
  (s3/delete-bucket client-config :bucket-name bucket))

(defn drop_bucket_nuclear
  [bucket]
  (while (not (empty? (list_bucket_contents bucket)))
    (do
      (doseq [object (list_bucket_contents bucket)]
        (drop_object bucket object))))
  (drop_bucket bucket)
  true)

(defn get_json
  ([bucket jsonpath]
   ;(log/infof "get_json request for bucket: %s  and file: %s" bucket filename)
   (try
     (let [s3object (s3/get-object client-config :bucket-name bucket :key (:name jsonpath) :prefix (:prefix jsonpath)) ; need to add a :prefix modeling ard storage on hsm
           s3content (clojure.java.io/reader (:object-content s3object))]
       (json/parse-stream s3content))
     (catch Exception e (log/errorf "Error retrieving data from object store: %s" e)
       false)))
  ([jsonpath]
   (get_json bucketname jsonpath)))

(defn get_url
  ([bucket filename expire]
   (.toString (s3/generate-presigned-url client-config bucket filename expire)))
  ([bucket filename]
   (let [today (jt/local-date)
         tomorrow (jt/plus today (jt/days 1))]
     (get_url bucket filename tomorrow)))
  ([filename]
   (get_url [bucketname filename])))

(defn init
  "Create bucket in object storage"
  []
  (create_bucket))

(mount/defstate bucket-init
  :start (init))
