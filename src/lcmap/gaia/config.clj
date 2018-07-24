(ns lcmap.gaia.config
  (:require [environ.core :as environ]))

(defn config
  []
  {:nemo_host (:nemo-host environ/env)  })
