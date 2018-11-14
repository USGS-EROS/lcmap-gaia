(ns lcmap.gaia.config
  (:require [environ.core :as environ]
            [cheshire.core :as json]))

(defn try-read
  [val]
  (try (read-string val)
       (catch Exception ex
         nil)))

(def lc_map (array-map
             :develop (or (try-read (:lc-develop environ/env)) 1)
             :ag      (or (try-read (:lc-ag      environ/env)) 2)
             :grass   (or (try-read (:lc-grass   environ/env)) 3)
             :tree    (or (try-read (:lc-tree    environ/env)) 4)
             :water   (or (try-read (:lc-water   environ/env)) 5)
             :wetland (or (try-read (:lc-wetland environ/env)) 6)
             :snow    (or (try-read (:lc-snow    environ/env)) 7)
             :barren  (or (try-read (:lc-barren  environ/env)) 8)
             :none    (or (try-read (:lc-none    environ/env)) 0)))

(def config
  {:nemo_host      (:nemo-host environ/env)
   :nemo_resource  (or (:nemo-resource environ/env) "/data")
   :fill_begin     (or (try-read (:fill-begin     environ/env)) true)
   :fill_end       (or (try-read (:fill-end       environ/env)) true)
   :fill_samelc    (or (try-read (:fill-samelc    environ/env)) true)
   :fill_difflc    (or (try-read (:fill-difflc    environ/env)) true)
   :fill_nodata    (or (try-read (:fill-nodata    environ/env)) true)
   :fill_nodataval (or (try-read (:fill-nodataval environ/env)) true)
   :lc_list        (vals lc_map)
   :lc_map         lc_map
   :lc_defaults (hash-map :lc_inbtw     (or (try-read (:lc-inbetween    environ/env)) 9)
                          :lc_insuff    (or (try-read (:lc-insufficient environ/env)) 10)
                          :lcc_growth   (or (try-read (:lc-growth       environ/env)) 151)
                          :lcc_decline  (or (try-read (:lc-decline      environ/env)) 152)
                          :lcc_nomodel  (or (try-read (:lc-nomodel      environ/env)) 201)
                          :lcc_forwards (or (try-read (:lc-forwards     environ/env)) 202)
                          :lcc_samelc   (or (try-read (:lc-samelc       environ/env)) 211)
                          :lcc_difflc   (or (try-read (:lc-difflc       environ/env)) 212)
                          :lcc_back     (or (try-read (:lc-back         environ/env)) 213)
                          :lcc_afterbr  (or (try-read (:lc-afterbreak   environ/env)) 214))})
