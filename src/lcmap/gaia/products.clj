(ns lcmap.gaia.products
  (:gen-class)
  (:require [clojure.tools.logging :as log]
            [clojure.stacktrace :as stacktrace]
            [clojure.string        :as string]
            [clojure.math.numeric-tower :as math]
            [clojure.walk          :refer [keywordize-keys]]
            [java-time             :as jt]
            [lcmap.gaia.file       :as file]
            [lcmap.gaia.util       :as util]
            [lcmap.gaia.config     :refer [config]]
            [lcmap.gaia.storage    :as storage]
            [lcmap.gaia.product-specs :as product-specs]
            [lcmap.gaia.nemo       :as nemo]))

(def product_abbreviations
  (hash-map "primary-landcover"              "LCPRI"
            "secondary-landcover"            "LCSEC"
            "primary-landcover-confidence"   "LCPCONF"
            "secondary-landcover-confidence" "LCSCONF"
            "annual-change"                  "LCACHG"
            "time-of-change"                 "SCTIME"
            "magnitude-of-change"            "SCMAG"
            "time-since-change"              "SCLAST"
            "magnitude-of-change"            "SCMQA"
            "length-of-segment"              "SCSTAB"))

(defn is-landcover
  [name]
  (let [abbr (get product_abbreviations name)]
    (try
      (some? (re-matches #"LC(.*)" abbr))
      (catch NullPointerException e
        false))))

(defn get-prefix
  [grid date tile]
  (let [hhh (subs tile 0 3)
        vvv (subs tile 3 6)
        year (first (string/split date #"-"))
        elements [year grid hhh vvv]]
    (string/join "/" elements)))

(defn map-path
  [tileid product date]
  ; LCMAP_<grid>_<6digit tileid>_<representative date>_<production date>_<CCDC version>_<product abbr>
  (let [grid      (:region config)
        repr_date (string/replace date "-" "")
        prod_date (string/replace (str (jt/local-date)) "-" "")
        ccd_ver   (:ccd_ver config)
        product_abbr (get product_abbreviations product)
        elements ["LCMAP" grid tileid repr_date prod_date ccd_ver product_abbr]
        name (str (string/join "-" elements) ".tif")
        prefix (get-prefix grid date tileid)
        url (storage/get_url storage/bucketname (str prefix "/" name))]
    {:name name :prefix prefix :url url}))

(defn ppath
  ([product x y tile date suffix]
   (let [grid      (:region config)
         name (->> [product x y date] (string/join "-") (#(str % suffix)))
         prefix (get-prefix grid date tile)]
     {:name name :prefix prefix}))
  ([product x y tile date]
   (ppath product x y tile date ".json")))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;    CHANGE PRODUCTS    ;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn product-exception-handler
  [exception product_name]
  (let [type_name (keyword (str product_name "-exception"))
        message (str "Error calculating " product_name)]
    (log/errorf (str message exception))
    (throw (ex-info message {:type type_name :cause :product-failure :exception exception}))))

(defn time-of-change
  "Return numeric day of year in which a break occurs"
  ([model query-day x y]
   (try
     (let [change-prob (:chprob model)
           break-day   (-> model (:bday) (util/to-javatime)) 
           query-year  (-> query-day (util/ordinal-to-javatime) (util/javatime-year))
           break-year  (-> break-day (util/javatime-year))
           response    #(hash-map :pixelx x :pixely y :val %)]
       (if (and (= query-year break-year) (= 1.0 change-prob))
         (-> break-day (util/javatime-day-of-year) response) 
         (-> 0 response)))
     (catch Exception e
       (product-exception-handler e "time-of-change"))))
  ([pixel_map pixel_models query-day]
   (let [values (map #(time-of-change % query-day (:px pixel_map) (:py pixel_map)) (:segments pixel_models))]
     (last (sort-by :val values)))))

(defn time-since-change
  "Return cumulative distance to previous break"
  ([model query-day x y]
   (try
     (let [change-prob (:chprob model)
           break-day   (-> model (:bday) (util/to-ordinal))
           day-diff    (- query-day break-day)
           distance    (if (and (= 1.0 change-prob) (pos? day-diff)) day-diff 0)]
       (hash-map :pixelx x :pixely y :val distance))
     (catch Exception e
       (product-exception-handler e "time-since-change"))))
  ([pixel_map pixel_models query-day]
   (let [values    (map #(time-since-change % query-day (:px pixel_map) (:py pixel_map)) (:segments pixel_models))
         non-zeros (filter (fn [i] (not (zero? (:val i)))) values)]
     (if (empty? non-zeros) 
       (first values)                    ; they are all 0, doesn't matter which 
       (first (sort-by :val non-zeros))) ; take the shortest distance
     )))

(defn magnitude-of-change
  "Return severity of spectral shift"
  ([model query-day x y]
   (try
     (let [change-prob (:chprob model)
           query-year  (-> query-day (util/ordinal-to-javatime) (util/javatime-year))
           break-year  (-> (:bday model) (util/to-javatime) (util/javatime-year))
           magnitudes  [(:grmag model) (:remag model) (:nimag model) (:s1mag model) (:s2mag model)]
           euc-norm    (math/sqrt (reduce + (map #(math/expt % 2) magnitudes)))
           response    #(hash-map :pixelx x :pixely y :val %)]
       (if (= true (= query-year break-year) (= 1.0 change-prob))
         (-> euc-norm (response))
         (-> 0 (response))))
     (catch Exception e
       (product-exception-handler e "magnitude-of-change"))))
  ([pixel_map pixel_models query-day]
   (let [values (map #(magnitude-of-change % query-day (:px pixel_map) (:py pixel_map)) (:segments pixel_models))]
     (last (sort-by :val values)))))

(defn length-of-segment
  "Return length of change segment in days"
  ([model query-day x y]
   (try
     (let [stability_begin (:stability_begin config)
           fill (- query-day (util/to-ordinal stability_begin))
           start-day (-> model (:sday) (util/to-ordinal))
           end-day   (-> model (:eday) (util/to-ordinal))
           diff      (if (> query-day end-day) (- query-day end-day) (- query-day start-day))
           value     (if (and (<= 0 diff) (< diff fill)) diff fill)]
       (hash-map :pixelx x :pixely y :val value))
     (catch Exception e
       (product-exception-handler e "length-of-segment"))))
  ([pixel_map pixel_models query-day]
   (let [values (map #(length-of-segment % query-day (:px pixel_map) (:py pixel_map)) (:segments pixel_models))]
     (first (sort-by :val values)))))

(defn curve-fit
  "Return Curve QA for point in time"
  ([model query-day x y]
   (try
     (let [curve-qa  (:curqa model)
           change-prob (:chprob model)
           query-year  (-> query-day (util/ordinal-to-javatime) (util/javatime-year))
           break-year  (-> (:bday model) (util/to-javatime) (util/javatime-year))
           value       (if (and (= change-prob 1.0) (= break-year query-year)) curve-qa 0)]
       (hash-map :pixelx x :pixely y :val value))
     (catch Exception e
       (product-exception-handler e "curve-fit"))))
  ([pixel_map pixel_models query-day]
   (let [values (map #(curve-fit % query-day (:px pixel_map) (:py pixel_map)) (:segments pixel_models))]
     (last (sort-by :val values)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;    CLASSIFICATION PRODUCTS    ;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn falls-between-eday-sday
  "Convenience function for returning pair of maps with true values for :follows_eday and :precedes_sday keys.
   Used with Reduce to identify maps in a list of sorted maps"
  [map_a map_b]
  (util/matching-keys map_a map_b :follows_eday :precedes_sday true))

(defn falls-between-bday-sday
  "Convenience function for returning pair of maps with true values for :follows_bday and :precedes_sday keys.
   Used with Reduce to identify maps in a list of sorted maps"
  [map_a map_b]
  (util/matching-keys map_a map_b :follows_bday :precedes_sday true))

(defn nbr
  "Return the Normalized Burn Ratio for a segment"
  [model]
  (let [sday   (-> model (:sday) (util/to-ordinal)) 
        eday   (-> model (:eday) (util/to-ordinal)) 
        niint  (:niint model)
        s1int  (:s1int model)
        nicoef (-> model (:nicoef) (first))
        s1coef (-> model (:s1coef) (first))
        nir_start  (+ niint (* sday nicoef))
        nir_end    (+ niint (* eday nicoef))
        swir_start (+ s1int (* sday s1coef))
        swir_end   (+ s1int (* eday s1coef))
        nbr_start  (float (/ (- nir_start swir_start) (+ nir_start swir_start)))
        nbr_end    (float (/ (- nir_end swir_end) (+ nir_end swir_end)))] 
    (- nbr_end nbr_start)))

(defn get-class 
 "Returns the class value given a collection of probabilities"
  ([probs rank]
   (let [sorted (reverse (sort probs)) 
         position (.indexOf probs (nth sorted rank))]
     (nth (:lc_list config) position)))
  ([probs]
   (get-class probs 0)))

(defn first-date-of-class
  "Returns the 'date' value from a collection of predictions for the first occurence of a given classification"
  [sorted_predictions class_val]
  (let [matching_predictions (filter (fn [i] (= class_val (get-class (:prob i)))) sorted_predictions)]
      (:date (first matching_predictions))))

(defn mean-probabilities
  "Returns a 1-d collection of mean probabilities given a collection of probabilities "
  [predictions]
  (let [probabilities (map :prob predictions)
        indexes (range 0 (count (first probabilities)))
        mean_map_fn (fn [i] (util/mean (map #(nth % i) probabilities)))]
    (map mean_map_fn indexes)))

(defn classify
  "Return the classification value for a single segment given a query_day and rank"
  ([model query_ord rank nbrdiff]
   (let [sorted_predictions (:probabilities model) ; sorted in characterize_segment
         first_class (-> sorted_predictions (first) (:prob) (get-class))
         last_class  (-> sorted_predictions (last)  (:prob) (get-class))
         grass_val (:grass (:lc_map config))
         tree_val  (:tree (:lc_map config))
         first_forest_date (-> sorted_predictions (first-date-of-class tree_val) (util/to-ordinal))   
         first_grass_date (-> sorted_predictions (first-date-of-class grass_val) (util/to-ordinal))    
         mean_probabilities (mean-probabilities sorted_predictions)]

     (cond
       ; nbr_slope > 0.05 and first_class is 'grass' and last is 'forest'
       (= true (> nbrdiff 0.05) (= grass_val first_class) (= tree_val last_class))
       (if (>= query_ord first_forest_date)
         (nth [tree_val grass_val] rank)
         (nth [grass_val tree_val] rank))

       ; nbr_slope < -0.05 and last class is grass and first class is forest
       (= true (< nbrdiff -0.05) (= tree_val first_class) (= grass_val last_class))
       (if (>= query_ord first_grass_date)
         (nth [grass_val tree_val] rank)
         (nth [tree_val grass_val] rank))
       :else ; calculate the mean across all probabilities for the segment, classify based on highest probability
       (get-class mean_probabilities rank))))
  ([model query_ord rank]
   (let [nbrdiff (nbr model)]
     (classify model query_ord rank nbrdiff))))

(defn characterize_segment
  "Return a hash-map characterizing details of the segment"
  [segment query_day probabilities rank]
  (let [sday (-> segment (:sday) (util/to-ordinal)) 
        eday (-> segment (:eday) (util/to-ordinal))
        bday (-> segment (:bday) (util/to-ordinal))
        nbrdiff (nbr segment)
        intersects        (<= sday query_day eday)
        precedes_sday     (< query_day sday)
        follows_eday      (> query_day eday)
        follows_bday      (>= query_day bday)
        between_eday_bday (<= eday query_day bday)
        growth  (> nbrdiff 0.05)
        decline (< nbrdiff -0.05)
        segment_probabilities (filter (fn [i] (= (-> i (:sday) (util/to-ordinal)) sday)) probabilities)
        sorted_probabilities  (util/sort-by-key segment_probabilities :date)
        classification (classify (merge segment {:probabilities sorted_probabilities}) query_day rank nbrdiff)]
    (hash-map :intersects     intersects
              :precedes_sday  precedes_sday
              :follows_eday   follows_eday
              :follows_bday   follows_bday
              :btw_eday_bday  between_eday_bday
              :sday           sday
              :eday           eday
              :bday           bday
              :growth         growth
              :decline        decline
              :probabilities  sorted_probabilities
              :classification classification)))

(defn landcover
  "Return the landcover value given the segments, probabilities, query_day and rank for a location"
  ([segments_probabilities query_day rank conf]
   (try
     (let [sorted_segments (util/sort-by-key (:segments segments_probabilities) :sday)
           probabilities   (:predictions segments_probabilities)
           characterized_segments (map #(characterize_segment % query_day probabilities rank) sorted_segments)
           first_start_day   (-> (first sorted_segments) (:sday) (util/to-ordinal))
           last_end_day      (-> (last sorted_segments)  (:eday) (util/to-ordinal))
           intersected_segment (first (filter :intersects characterized_segments))
           eday_bday_model   (first (filter :btw_eday_bday characterized_segments))
           between_eday_sday (reduce falls-between-eday-sday characterized_segments)
           between_bday_sday (reduce falls-between-bday-sday characterized_segments)]

       (cond
         ; query date precedes first segment start date and fill_begin config is true
         (= true (< query_day first_start_day) (:fill_begin conf))
         (:classification (first characterized_segments)) ; return value of the first segment

         ; query date precedes first segment start date
         (= true (< query_day first_start_day))
         (:lc_insuff (:lc_defaults conf)) ; return lc_insuff value from lc_defaults config

         ; query date follows last segment end date and fill_end config is true
         (= true (> query_day last_end_day) (:fill_end conf))
         (:classification (last characterized_segments)) ; return value of the last segment
         
         ; query date follows last segment end date
         (= true (> query_day last_end_day))
         (:lc_insuff (:lc_defaults conf)) ; return the lc_insuff value from the lc_defaults config

         ; query date falls between a segments start date and end date
         (not (nil? intersected_segment))
         (:classification intersected_segment) ; return the class value for the intercepted model

         ; query date falls between segments of same landcover classification and fill_samelc config is true
         (= true (:fill_samelc conf) (= (:classification (first between_eday_sday)) (:classification (last between_eday_sday))))
         (:classification (last between_eday_sday)) ; return the value from the last model from the pair of models the query date fell between

         ; query date falls between one segments break date and the following segments start date and fill_difflc config is true
         (= true (:fill_difflc conf) (not (map? between_bday_sday)))
         (:classification (last between_bday_sday )) ; return the value from the last model from the pair of models the query date fell between

         ; query date falls between a segments end date and break date and fill_difflc config is true
         (= true (:fill_difflc conf) (not (nil? eday_bday_model)))
         (:classification eday_bday_model) ; return the value from the model where the query date intersected the end date and break date

         :else ; finally as a last resort return the lc_inbtw value from the configuration
         (:lc_inbtw conf)))

     (catch Exception e
       (product-exception-handler e "landcover"))))
  ([segments_probabilities query_day rank] ; enable passing in the configuration
   (landcover segments_probabilities query_day rank config)))

(defn primary-landcover
  "Return the highest landcover class value"
  [pixel_coords segments_probabilities query_day]
  (let [value (landcover segments_probabilities query_day 0)]
    (hash-map :pixelx (:px pixel_coords) :pixely (:py pixel_coords) :val value)))

(defn secondary-landcover
  "Return the second highest landcover class value"
  [pixel_coords segments_probabilities query_day]
  (let [value (landcover segments_probabilities query_day 1)]
    (hash-map :pixelx (:px pixel_coords) :pixely (:py pixel_coords) :val value)))

(defn annual-change
  "Return the change in landcover from the provided year, to the previous year"
  [pixel_coords segments_probabilities query_ord]
  (let [previous_query_day (util/subtract_year query_ord)
        previous_value (landcover segments_probabilities previous_query_day 0)
        latest_value (landcover segments_probabilities query_ord 0)
        response_template (hash-map :pixelx (:px pixel_coords) :pixely (:py pixel_coords))]

    (if (= previous_value latest_value)
      (merge response_template {:val latest_value})
      (merge response_template {:val (util/concat_ints previous_value latest_value)}))))

(defn confidence
  "Return the landcover confidence value given the segments, probabilities, query_day and rank for a location"
  [segments_probabilities query_day rank]
  (try
    (let [sorted_segments (util/sort-by-key (:segments segments_probabilities) :sday)
          probabilities   (:predictions segments_probabilities)
          characterized_segments (map #(characterize_segment % query_day probabilities rank) sorted_segments)
          first_start_day   (-> (first sorted_segments) (:sday) (util/to-ordinal))
          last_end_day      (-> (last sorted_segments)  (:eday) (util/to-ordinal))
          intersected_segment (first (filter :intersects characterized_segments))
          eday_bday_model   (first (filter :btw_eday_bday characterized_segments))
          between_eday_sday (reduce falls-between-eday-sday characterized_segments)
          between_bday_sday (reduce falls-between-bday-sday characterized_segments)]

      (cond
        ; query date preceds first segment start date
        (= true (< query_day first_start_day))
        (:lcc_back (:lc_defaults config)) ; return lcc_back value from lc_defaults config

        ; query date follows last segment end date and change prob == 1
        (= true (> query_day last_end_day) (= 1 (int (:chprob (last sorted_segments)))))
        (:lcc_afterbr (:lc_defaults config)) ; return the lcc_afterbr value from the lc_defaults config

        ; query date follows last segment end date
        (= true (> query_day last_end_day))
        (:lcc_forwards (:lc_defaults config)) ; return the lcc_forwards value from the lc_defaults config

        ; query date falls between a segments start date and end date and growth is true
        (= true (not (nil? intersected_segment)) (:growth intersected_segment))
        (:lcc_growth (:lc_defaults config)) ; return lcc_growth value from lc_defaults config

        ; query date falls between a segments start date and end date and decline is true
        (= true (not (nil? intersected_segment)) (:decline intersected_segment))
        (:lcc_decline (:lc_defaults config)) ; return lcc_decline value from lc_defaults config

        ; query date falls between a segments start date and end date
        (not (nil? intersected_segment))
        (util/scale-value (nth (:prob (last (:probabilities intersected_segment))) rank))

        ; query date falls between segments of same landcover classification
        (= true (= (:classification (first between_eday_sday)) (:classification (last between_eday_sday))))
        (:lcc_samelc (:lc_defaults config)) ; return lcc_samelc from lc_defaults config

        ; query date falls between segments with different landcover classifications
        (= 2 (count between_eday_sday))
        (:lcc_difflc (:lc_defaults config)) ; return lcc_difflc from lc_defaults config

        :else ; mapify returns ValueError
        (:none (:lc_map config))))

    (catch Exception e
      (product-exception-handler e "confidence"))))

(defn primary-landcover-confidence
  "Return the landcover probability for the highest landcover class value"
  [pixel_coords segments_probabilities query_day]
  (let [value (confidence segments_probabilities query_day 0)]
    (hash-map :pixelx (:px pixel_coords) :pixely (:py pixel_coords) :val value)))

(defn secondary-landcover-confidence
  "Return the landcover probability for the 2nd highest landcover class value"
  [pixel_coords segments_probabilities query_day]
  (let [value (confidence segments_probabilities query_day 1)]
    (hash-map :pixelx (:px pixel_coords) :pixely (:py pixel_coords) :val value)))

(defn pixel_map
  "Return hash-map keyed by pixelx and pixely with a hash-map value for :segments and :predictions"
  [inputs]
  (let [pixelx      (first (:pixelxy inputs))
        pixely      (last  (:pixelxy inputs))
        segments    (get (:segments inputs)    (:pixelxy inputs)) 
        predictions (get (:predictions inputs) (:pixelxy inputs))]
    (hash-map {:px pixelx :py pixely} (hash-map :segments segments :predictions predictions))))

(defn pixel_groups
  [injson]
  (let [juxt_fn (util/variable-juxt [:px :py])]
    (group-by juxt_fn injson)))

(defn flatten_product_data
  "Return a flat list of product values given a collection of hash-maps
  for every pixel in a chip, [{:pixely 3159045, :pixelx -2114775, :val 6290},...] ...]"
  [product_value_collection]
  (let [; group product coll by row
        row_groups (util/coll-groups product_value_collection [:pixely]) 
        ; sort row group values by pixelx ascending 
        sort-pixelx-fn (fn [i] (hash-map (:pixely (first i)) (sort-by :pixelx (last i))))
        sorted-x-vals (map sort-pixelx-fn row_groups)
        ; sort the rows by the pixely key ascending
        sorted-y-rows (sort-by (fn [i] (first (keys i))) > sorted-x-vals)
        ; finally, flatten to a one dimensional list
        flattened (util/flatten-vals sorted-y-rows :val)]
    flattened))

(defn data
  "Returns a 1-d collection of product values"
  [segments_json predictions_json product_type queryday]
  (let [segments    (-> segments_json (keywordize-keys) (product-specs/segment_coll_check) (pixel_groups))
        predictions (or (some-> predictions_json (not-empty) (keywordize-keys) (product-specs/prediction_coll_check) (pixel_groups)) [])
        product_fn  (->> product_type (product-specs/product_type_check) (str "lcmap.gaia.products/") (symbol) (resolve))
        query_ord   (-> queryday (product-specs/date_fmt_check) (util/to-ordinal))
        per_pixel_inputs (map #(pixel_map {:pixelxy % :segments segments :predictions predictions}) (keys segments))
        per_pixel_values (map #(product_fn (first (keys %)) (first (vals %)) query_ord) per_pixel_inputs)]
    (-> per_pixel_values (flatten_product_data) (product-specs/output_check))))

(defn persist
  [product cx cy tile query_day segments predictions]
  (try
    (let [values (data segments predictions product query_day) 
          out_path (ppath product cx cy tile query_day)
          out_data {"x" cx "y" cy "values" values}]
      (log/infof "storing : %s" (:name out_path))
      (storage/put_json out_path out_data)
      {:cx cx :cy cy :date query_day :status "success"})

    (catch Exception e 
      (log/errorf "Exception in persist - cx: %s  cy: %s  product: %s date: %s exception-message: %s exception-data: %s" 
                  cx cy product query_day (.getMessage e) (ex-data e))
           {:cx cx :cy cy :date query_day :product product :status "fail" :message (str (.getMessage e) " - " (ex-data e))})))

(defn generation
  [{dates :dates cx :cx cy :cy product :product tile :tile :as all}]
  (try
    (let [segments (nemo/segments cx cy)
          predictions (if (is-landcover product) (nemo/predictions cx cy) []) ; predictions are not required for change products. don't make unnecessary http requests
          persist #(persist product cx cy tile % segments predictions)
          results (pmap persist dates)
          failures (->> results 
                        (filter (fn [i] (= "fail" (:status i)))) 
                        (map (fn [i] {(:date i) (:message i)})))]
      {:failures failures :product product :cx cx :cy cy :dates dates})
    (catch Exception e
      (log/errorf "Exception in products/generation ! args: %s -  message: %s - data: %s - stacktrace: %s" all (.getMessage e) (ex-data e) (-> e stacktrace/print-stack-trace with-out-str))
      (throw (ex-info "Exception in product generation" {:data (ex-data e) :error-message (.getMessage e) :args all})))))
