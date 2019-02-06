(ns lcmap.gaia.storage
  (:require [lcmap.gaia.config :refer [config]]
            [amazonica.aws.s3 :as s3]))

(def client-config
  {:access-key (:s3-access-key config)
   :secret-key (:s3-secret-key config)
   :endpoint   (:s3-endpoint config)
   :client-config {:path-style-access-enabled true}})

(defn list_buckets
  []
  (s3/list-buckets client-config))

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


;; s3test.core=> (def x (assoc cred :client-config {:path-style-access-enabled true}))
;; #'s3test.core/x
;; s3test.core=> (s3/list-buckets x)
;; [{:creation-date #object[org.joda.time.DateTime 0x224b55a0 "2018-12-17T08:06:24.547-06:00"], :owner {:display-name "cgaustin", :id "cgaustin"}, :name "clayfoo"}]
;; s3test.core=> 


;; s3test.core=> (s3/create-bucket x "jack")
;; {:name "jack"}
;; s3test.core=> (s3/list-buckets x)


;; s3test.core=> (def upload-file (java.io.File. "upload.txt"))
;; #'s3test.core/upload-file
;; s3test.core=> (spit upload-file "hello world")
;; nil
;; s3test.core=> (s3/put-object x :bucket-name "jack" :key "jenny" :file upload-file)
