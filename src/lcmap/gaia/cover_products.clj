(ns lcmap.gaia.cover-products
  (:gen-class)
  (:require [clojure.core.async    :as async]
            [clojure.math.numeric-tower :as math]
            [clojure.math.combinatorics :as combo]
            [clojure.string        :as string]
            [clojure.tools.logging :as log]
            [java-time             :as jt]
            [lcmap.gaia.config     :refer [config]]
            [lcmap.gaia.file       :as file]
            [lcmap.gaia.gdal       :as gdal]
            [lcmap.gaia.product-specs :as product-specs]
            [lcmap.gaia.storage    :as storage]
            [lcmap.gaia.util       :as util]
            [cheshire.core         :as json]))

(defn product-exception-handler
  [exception product_name pixel]
  (let [msg (format "problem calculating %s with pixel %s: %s" product_name pixel (.getMessage exception))]
    (log/error msg)
    (throw (ex-info msg {:type "data-generation-error" :message msg} (.getCause exception)))))

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

(defn normalized-burn-ratio
  "Return the Normalized Burn Ratio for a segment"
  [model sday eday]
  (let [niint  (get model "niint")
        s1int  (get model "s1int")
        nicoef (first (get model "nicoef"))
        s1coef (first (get model "s1coef"))
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
   (try
     (let [sorted (reverse (sort probs)) 
           position (.indexOf probs (nth sorted rank))]
       (nth (:lc_list config) position))
     (catch IndexOutOfBoundsException e ; the probs collection is empty
       (:none (:lc_map config)))))      ; return configured value for None
  ([probs]
   (get-class probs 0)))

(defn first-date-of-class
  "Returns the 'date' value from a collection of predictions for the first occurence of a given classification"
  [sorted_predictions class_val]
  (let [matching_predictions (filter (fn [i] (= class_val (get-class (get i "prob")))) sorted_predictions)]
      (get (first matching_predictions) "pday")))

(defn mean-probabilities
  "Returns a 1-d collection of mean probabilities given a collection of probabilities "
  [predictions]
  (let [probabilities (map #(get % "prob") predictions)
        indexes (range 0 (count (first probabilities)))
        mean_fn (fn [i] (util/mean (map #(nth % i) probabilities)))]
    (map mean_fn indexes)))

(defn scaled-probabilities 
  [probabilities rank] 
  (let [operation (comp util/scale-value #(nth % rank) reverse sort mean-probabilities)] 
    (operation probabilities)))

(defn class-details
  [predictions query_date rank burn_ratio]
  (let [first_class (get-class (get (first predictions) "prob"))
        last_class  (get-class (get (last predictions) "prob"))
        grass (:grass (:lc_map config))
        tree  (:tree  (:lc_map config))
        first_forest_date (util/to-ordinal (first-date-of-class predictions tree))   
        first_grass_date  (util/to-ordinal (first-date-of-class predictions grass))    
        growth (and (> burn_ratio 0.05) (= grass first_class) (= tree last_class))
        decline (and (< burn_ratio -0.05) (= tree first_class) (= grass last_class))
        probabilities (mean-probabilities predictions)
        class (get-class probabilities rank)]
    (hash-map
     :first_class first_class
     :last_class  last_class
     :first_forest_date first_forest_date
     :first_grass_date first_grass_date
     :mean_probabilities probabilities
     :growth growth
     :decline decline
     :class class)))

(defn classify
  "Return the classification value for a single segment given a query_day and rank"
  [predictions query_date rank burn_ratio]
  (let [grass (:grass (:lc_map config))
        tree  (:tree  (:lc_map config))
        details (class-details predictions query_date rank burn_ratio)]

    (cond
     (:growth details)
     (if (>= query_date (:first_forest_date details))
       (nth [tree grass] rank)
       (nth [grass tree] rank))

     (:decline details)
     (if (>= query_date (:first_grass_date details))
       (nth [grass tree] rank)
       (nth [tree grass] rank))

     :else
     (:class details))))

(defn characterize-segment
  "Return a hash-map characterizing details of the segment"
  [segment query_day probabilities]
  (let [sday ((comp util/to-ordinal #(get % "sday")) segment)
        eday ((comp util/to-ordinal #(get % "eday")) segment)
        bday ((comp util/to-ordinal #(get % "bday")) segment)
        chprob (get segment "chprob")
        burn_ratio (normalized-burn-ratio segment sday eday)
        intersects        (<= sday query_day eday)
        precedes_sday     (< query_day sday)
        follows_eday      (> query_day eday)
        follows_bday      (>= query_day bday)
        between_eday_bday (<= eday query_day bday)
        ordinal_sday #(util/to-ordinal (get % "sday"))
        probability_reducer (fn [coll p] (if (= (ordinal_sday p) sday) (conj coll p) coll)) ; if prediction sday == segment sday, keep it
        segment_probabilities (reduce probability_reducer [] probabilities)
        sorted_probabilities  (util/sort-by-key segment_probabilities "pday")
        primary_classification   (classify sorted_probabilities query_day 0 burn_ratio)
        secondary_classification (classify sorted_probabilities query_day 1 burn_ratio)
        class_details (class-details sorted_probabilities query_day 0 burn_ratio)]
    (hash-map :intersects      intersects
              :precedes_sday   precedes_sday
              :follows_eday    follows_eday
              :follows_bday    follows_bday
              :btw_eday_bday   between_eday_bday
              :sday            sday
              :eday            eday
              :bday            bday
              :growth          (:growth class_details)
              :decline         (:decline class_details)
              :chprob          chprob
              :probabilities   sorted_probabilities
              :primary_class   primary_classification
              :secondary_class secondary_classification)))

(defn same-class
  [segments class_key]
  (let [seg_one (first segments)
        seg_two (last segments)] 
    (= (class_key seg_one) (class_key seg_two))))

(defn landcover
  "Return the landcover value given the segments, probabilities, query_day and rank for a location"
  ([segments query_date conf]
   (try
     (let [first_start_day     (:sday (first segments))
           last_end_day        (:eday (last segments))
           intersected_segment (first (filter :intersects segments))
           eday_bday_model     (first (filter :btw_eday_bday segments))
           between_eday_sday   (reduce falls-between-eday-sday segments)
           between_bday_sday   (reduce falls-between-bday-sday segments)
           ;class_key           (if (= 0 rank) :primary_class :secondary_class)
           key_response        #(hash-map :primary (:primary_class %) :secondary (:secondary_class %))]

       (cond
         ; query date precedes first segment start date
         (< query_date first_start_day)
         (key_response (first segments))  ; (class_key (first segments)) ; return value of the first segment

         ; query date follows last segment end date
         (> query_date last_end_day)
         (key_response (last segments)) ;(class_key (last segments)) ; return value of the last segment
         
         ; query date falls between a segments start date and end date
         (not (nil? intersected_segment))
         (key_response intersected_segment) ; (class_key intersected_segment) ; return the class value for the intercepted model

         ;; ; query date falls between segments of same landcover classification and fill_samelc config is true
          (and (:fill_samelc conf) (= (class_key (first between_eday_sday)) (class_key (last between_eday_sday))))
          (key_response (last between_eday_sday)) ;(class_key (last between_eday_sday)) ; return the value from the last model from the pair of models the query date fell between

          ; both are the same

          ; either primary or secondary is the same
          ; for the primary/secondary thats the same, take the primary/secondary of the (last between_eday_sday)
          ; for the other: (if (:fill_difflc conf)
          ;                  


          ; neither are the same ..

         ; query date falls between one segments break date and the following segments start date and fill_difflc config is true
         (and (:fill_difflc conf) (not (map? between_bday_sday)))
         (key_response (last between_bday_sday)) ;(class_key (last between_bday_sday)) ; return the value from the last model from the pair of models the query date fell between

         ; query date falls between a segments end date and break date and fill_difflc config is true
         (and (:fill_difflc conf) (not (nil? eday_bday_model)))
         (key_response eday_bday_model) ;(class_key eday_bday_model) ; return the value from the model where the query date intersected the end date and break date

         :else ; we need to throw an exception
         (throw (ex-info (format "Landcover value calculation problem, unclassifiable pixel %s" (:pixelxy segments)) 
                         {:type "data-generation-error"}))))

     (catch Exception e
       (product-exception-handler e "landcover" segments))))
  ([characterized_pixel date] ; enable passing in the configuration
   (landcover characterized_pixel date config)))

(defn change
  "Return the change in landcover from the provided year, to the previous year"
  [current_landcover previous_landcover]
  (if (= current_landcover previous_landcover)
    current_landcover
    (util/concat_ints previous_landcover current_landcover)))


(defn confidence
  "Return the landcover confidence value given the segments, probabilities, query_day and rank for a location"
  ([segments query_date conf]
   (try
     (let [[px py]             (:pixelxy segments)
           first_start_day     (:sday (first segments))
           last_end_day        (:eday (last segments))
           intersected_segment (first (filter :intersects segments))
           eday_bday_model     (first (filter :btw_eday_bday segments))
           between_eday_sday   (reduce falls-between-eday-sday segments)
           between_bday_sday   (reduce falls-between-bday-sday segments)
           key_response        #(hash-map :primary % :secondary %)
           ;class_key           (if (= 0 rank) :primary_class :secondary_class)
           ]

       (cond
        ; first segment wasn't classifiable, and query date falls between segments
        (and (= 2 (count between_eday_sday)) (empty? (:probabilities (first between_eday_sday))))
        (key_response (:lcc_back (:lc_defaults conf)))     ;(:lcc_back (:lc_defaults conf))

        ; last segment wasn't classifiable, and query date falls between segments
        (and (= 2 (count between_eday_sday)) (empty? (:probabilities (last between_eday_sday))))
        (key_response (:lcc_afterbr (:lc_defaults conf)))  ;(:lcc_afterbr (:lc_defaults conf))

        ; query date precedes first segment start date
        (< query_date first_start_day)
        (key_response (:lcc_back (:lc_defaults conf)))  ;(:lcc_back (:lc_defaults conf)) ; return lcc_back value from lc_defaults config

        ; query date follows last segment end date and change prob == 1
        (and (> query_date last_end_day) (= 1 (int (:chprob (last segments)))))
        (key_response (:lcc_afterbr (:lc_defaults conf)))  ;(:lcc_afterbr (:lc_defaults conf)) ; return the lcc_afterbr value from the lc_defaults config

        ; query date follows last segment end date and last segment has no predictions
        (and (> query_date last_end_day) (empty? (:probabilities (last segments))))
        (key_response (:lcc_afterbr (:lc_defaults conf)))  ;(:lcc_afterbr (:lc_defaults conf)) ; return the lcc_afterbr value from the lc_defaults config

        ; query date follows last segment end date
        (> query_date last_end_day)
        (key_response (:lcc_forwards (:lc_defaults conf)))  ;(:lcc_forwards (:lc_defaults conf)) ; return the lcc_forwards value from the lc_defaults config

        ; query date falls between a segments start date and end date and growth is true
        (and (not (nil? intersected_segment)) (:growth intersected_segment))
        (key_response (:lcc_growth (:lc_defaults conf)))  ;(:lcc_growth (:lc_defaults conf)) ; return lcc_growth value from lc_defaults config

        ; query date falls between a segments start date and end date and decline is true
        (and (not (nil? intersected_segment)) (:decline intersected_segment))
        (key_response (:lcc_decline (:lc_defaults conf)))  ;(:lcc_decline (:lc_defaults conf)) ; return lcc_decline value from lc_defaults config

        ; query date falls between a segments start date and end date
        (and (not (nil? intersected_segment)) (not (empty? (:probabilities intersected_segment))))
        (hash-map :primary (scaled-probabilities (:probabilities intersected_segment) 0) 
                  :secondary (scaled-probabilities (:probabilities intersected_segment) 1))
        ;(-> (:probabilities intersected_segment) mean-probabilities sort reverse (nth rank) util/scale-value)

        (= 2 (count between_eday_sday))
        (hash-map :primary (if (same-class between_eday_sday :primary_class) 
                             (:lcc_samelc (:lc_defaults conf))   ; query date falls between segments of same landcover classification
                             (:lcc_difflc (:lc_defaults conf)))  ; query date falls between segments with different landcover classifications
                  :secondary (if (same-class between_eday_sday :secondary_class) 
                             (:lcc_samelc (:lc_defaults conf))   ; query date falls between segments of same landcover classification
                             (:lcc_difflc (:lc_defaults conf)))) ; query date falls between segments with different landcover classifications

        ;; ; query date falls between segments of same landcover classification
        ;; (= (class_key (first between_eday_sday)) (class_key (last between_eday_sday)))
        ;; (:lcc_samelc (:lc_defaults conf)) ; return lcc_samelc from lc_defaults config

        ;; ; query date falls between segments with different landcover classifications
        ;; (= 2 (count between_eday_sday))
        ;; (:lcc_difflc (:lc_defaults conf)) ; return lcc_difflc from lc_defaults config

        :else ; we need to throw and exception
        (throw (ex-info (format "Confidence value calculation problem, pixel %s" (:pixelxy segments)) 
                         {:type "data-generation-error"}))))
     (catch Exception e
       (product-exception-handler e "confidence" segments))))
  ([characterized_pixel date]
   (confidence characterized_pixel date config)))

(defn characterize-inputs
  "Return a hash-map characterizing details of the segment"
  [pixelxy inputs query_day]
  (let [segments_valid    (product-specs/segments-valid? (:segments inputs))
        predictions_valid (product-specs/predictions-valid? (:predictions inputs))
        previous_day      (util/subtract_year query_day)
        characterize_fn   (fn [d] #(characterize-segment % d (:predictions inputs)))
        kw_query_day      (keyword (str query_day))
        kw_previous_day   (keyword (str previous_day))]
    (if (and segments_valid predictions_valid)
      (hash-map :pixelxy pixelxy
                :current {kw_query_day (map (characterize_fn query_day) (:segments inputs))}
                :previous {kw_previous_day (map (characterize_fn previous_day) (:segments inputs))})
      (hash-map :pixelxy pixelxy 
                :current {kw_query_day []} 
                :previous {kw_previous_day []}))))

(defn products 
  [characterized_pixel]
  (let [kw_date            (first (keys (:current characterized_pixel)))
        kw_previous        (first (keys (:previous characterized_pixel)))
        date               (read-string (name kw_date)) 
        previous_date      (read-string (name kw_previous)) 
        current_character  (get-in characterized_pixel [:current kw_date])
        previous_character (get-in characterized_pixel [:previous kw_previous]) 
        [px py] (:pixelxy characterized_pixel)
        none    (:none (:lc_map config))
        nomodel (:lcc_nomodel (:lc_defaults config))
        good_data (not (empty? current_character))
        segments (filter (fn [i] (not (empty? (:probabilities i)))) current_character)
        previous_segments (filter (fn [i] (not (empty? (:probabilities i)))) previous_character)

        ;primary_landcover    (if good_data (landcover  segments date 0) none)
        ;secondary_landcover  (if good_data (landcover  segments date 1) none)
        lc (landcover segments date)
        primary_landcover   (if good_data (:primary   lc) none)
        secondary_landcover (if good_data (:secondary lc) none)

        ;primary_confidence   (if good_data (confidence segments date 0) nomodel) 
        ;secondary_confidence (if good_data (confidence segments date 1) nomodel)
        lc_conf (confidence segments date)
        primary_confidence   (if good_data (:primary   lc_conf) none)
        secondary_confidence (if good_data (:secondary lc_conf) none)

        previous_landcover   (if (not (empty? previous_character)) (landcover previous_segments previous_date) none)
        annual_change        (if good_data (change primary_landcover (:primary previous_landcover)) none)]

    (hash-map :px px :py py :date date
              :values {:primary-landcover primary_landcover 
                       :secondary-landcover secondary_landcover
                       :primary-confidence primary_confidence
                       :secondary-confidence secondary_confidence
                       :annual-change annual_change})))

(defn start-consumers
  [in-chan out-chan operation]
  (let [chunk_size  (:product_instance_count config)]
    (dotimes [_ chunk_size]
      (async/thread
        (while true
          (let [input (async/<!! in-chan)
                result (operation input)]
            (async/>!! out-chan result)))))))

(defn generate
  [{dates :dates cx :cx cy :cy tile :tile :as all}]
  (try
    (let [grouped_segments     (storage/grouped-segments cx cy)
          grouped_predictions  (storage/grouped-predictions cx cy)
          in-chan              (async/chan)
          out-chan             (async/chan)
          ; pixel coordinates
          pixels               (keys grouped_segments) 
          pixel_hash-map       #(hash-map % {:segments (get grouped_segments %) :predictions (get grouped_predictions %)})
          ; segments and preds grouped by px py
          pixel_inputs         (into {} (map pixel_hash-map pixels)) 
          ordinal_dates        (doall (map util/to-ordinal dates))
          ; collection of all date and pixel coordinate combinations 
          dates_pixels         (combo/cartesian-product ordinal_dates pixels)
          ; characterize the segments, then calculate product values
          operation            (comp products #(characterize-inputs (last %) (get pixel_inputs (last %)) (first %)))
          ; items put on in-chan, will be operated on by the operation fn, results put onto the out-chan
          consumers            (start-consumers in-chan out-chan operation)
          output_fn            (fn [i] (let [result (async/<!! out-chan)] result))]

      (async/go
        (doseq [dp dates_pixels] ; dp looks like (726283 [1784625 2414775])
          (async/>!! in-chan dp)))

      (let [results (doall (map output_fn dates_pixels))
            date_groups (group-by :date results)]
        (doseq [dg date_groups
                :let [date (util/to-yyyy-mm-dd (first dg))
                      path (storage/ppath "cover" cx cy tile date)
                      values (util/flatten-values (second dg))]]
          ;(util/with-retry (storage/put_json path values))
          ))

      {:products "cover" :cx cx :cy cy :dates dates})
    (catch Exception e
      (let [msg (format "problem generating cover products for %s: %s" all (.getMessage e))]
        (log/error msg)
        (throw (ex-info msg {:type "data-generation-error" :message msg} (.getCause e)))))))
