(ns lcmap.gaia.config
  (:require [environ.core :as environ]))

(def config
  {:nemo_host (:nemo-host environ/env)
   :nemo_resource (or (:nemo-resource environ/env) "/data")})
