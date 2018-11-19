(defproject faux "0.1.0-SNAPSHOT"
  :description "Dumb JSON file server"
  :url "fake-gnus"
  :license {:name "Unlicense"
            :url ""}
  :dependencies [[org.clojure/clojure       "1.9.0"]
                 [cheshire                  "5.8.0"]
                 [compojure                 "1.6.1"]
                 [com.cemerick/pomegranate  "1.0.0" :exclusions [[org.slf4j/jcl-over-slf4j]
                                                                 [org.apache.httpcomponents/httpcore]
                                                                 [org.slf4j/slf4j-api]]]
                 [environ                   "1.1.0"]
                 [http-kit                  "2.2.0"]
                 [http-kit.fake             "0.2.1"]
                 [ring                      "1.6.3"]
                 [ring/ring-defaults        "0.3.1"]
                 [ring/ring-json            "0.4.0"]
                 [ring/ring-jetty-adapter   "1.6.3"]
                 [ring/ring-mock            "0.3.2"]]

  :plugins [[lein-environ "1.1.0"]]
  :profiles {:dev     {:dependencies [[org.clojure/test.check "0.9.0"]]}
             :uberjar {:omit-source true
                       :aot :all}}
  :main faux.main)
