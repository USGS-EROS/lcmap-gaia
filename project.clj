(defproject gaia "0.1.0-SNAPSHOT"
  :description "Level 2 CCDC Product Generation"
  :url "http://github.com/usgs-eros/lcmap-gaia"
  :license {:name "Unlicense"
            :url ""}
  :dependencies [[cheshire                  "5.8.0"]
                 [clojure.java-time         "0.3.2"]
                 [environ                   "1.1.0"]
                 [org.clojure/clojure       "1.9.0"]
                 [org.clojure/core.async    "0.3.443"]
                 [org.clojure/tools.logging "0.4.0"]]
  :plugins [[lein-environ "1.1.0"]]
  :profiles {:dev     {:resource-paths ["dev"]}
             :repl    {:resource-paths ["dev"]
		       :dependencies [[cider/cider-nrepl "0.15.1"]]}
             :test    {:resource-paths ["test" "test/resources"]}
             :uberjar {:omit-source true
                       :aot :all}}
  :main lcmap.gaia.main)
