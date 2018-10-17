(ns lcmap.gaia.test-resources
  (:require [lcmap.gaia.file :as file]
            [lcmap.gaia.util :as util]))

(def chip_data (file/read-json "resources/cx-2115585_cy3119805_segment.json"))
(def pixel_segments (util/pixel-groups chip_data))
(def grid_data (file/read-json "resources/grid.json"))
