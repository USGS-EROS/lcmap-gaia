(ns lcmap.gaia.test-resources
  (:require [lcmap.gaia.file :as file]
            [lcmap.gaia.util :as util]))

(def chipdata (file/read-json "resources/y3161805_x-2115585_nodates.json"))
(def pixel_segments (util/pixel-groups chipdata))
