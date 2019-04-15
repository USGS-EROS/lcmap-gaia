(defproject gaia "0.1.3"
  :description "Level 2 CCDC Product Generation"
  :url "http://github.com/usgs-eros/lcmap-gaia"
  :license {:name "Unlicense"
            :url ""}
  :dependencies [[org.clojure/clojure       "1.9.0"]
                 [org.clojure/core.async    "0.3.443"]
                 [org.clojure/math.numeric-tower "0.0.4"]
                 [org.clojure/tools.logging "0.4.0"]
                 [org.slf4j/slf4j-log4j12 "1.7.21"]
                 [log4j/log4j "1.2.17" :exclusions [javax.mail/mail
                                                    javax.jms/jms
                                                    com.sun.jmdk/jmxtools
                                                    com.sun.jmx/jmxri]]
                 [cheshire                  "5.8.0"]
                 [clojure.java-time         "0.3.2"]
                 [compojure                 "1.6.1"]
                 [environ                   "1.1.0"]
                 [http-kit                  "2.2.0"]
                 [http-kit.fake             "0.2.1"]
                 [ring                      "1.6.3"]
                 [ring/ring-defaults        "0.3.1"]
                 [ring/ring-json            "0.4.0"]
                 [ring/ring-jetty-adapter   "1.6.3"]
                 [ring/ring-mock            "0.3.2"]
                 [org.gdal/gdal             "2.2.0"]
                 [mount                     "0.1.12"]
                 [amazonica                 "0.3.139"]]

  :plugins [[lein-environ "1.1.0"]]
  :profiles {:dev     {:dependencies [[org.clojure/test.check "0.9.0"]]}
             :test    {:resource-paths ["test" "test/resources"] :dependencies [[org.clojure/test.check "0.9.0"]]}
             :uberjar {:omit-source true
                       :aot :all}}
  :main lcmap.gaia.main)
