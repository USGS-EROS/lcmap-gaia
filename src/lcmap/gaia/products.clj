(ns lcmap.gaia.products
  (:gen-class)
  (:require [clojure.tools.logging :as log]
            [clojure.string        :as string]
            [clojure.math.numeric-tower :as math]
            [lcmap.gaia.file       :as file]
            [lcmap.gaia.util       :as util]
            [lcmap.gaia.config     :refer [config]]))

(defn sort-by-sday [models] (sort-by (fn [i] (get i "sday")) models))

(defn product-name
  [model product fmt]
  (let [chipx (get model "cx")
        chipy (get model "cy")
        name  (string/join "_" [product chipx chipy])]
    (str name "." fmt)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;    CHANGE PRODUCTS    ;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn time-of-change
  "Return numeric day of year in which a break occurs"
  ([model query-day x y]
   (let [change-prob (get model "chprob")
         break-day   (-> (get model "bday") (util/to-ordinal)) 
         query-year  (-> query-day (util/to-javatime) (util/javatime-year))
         break-year  (-> break-day (util/ordinal-to-javatime) (util/javatime-year))
         response    #(hash-map :pixelx x :pixely y :val %)]
     (if (= true (= query-year break-year) (= 1.0 change-prob))
       (-> break-day (util/ordinal-to-javatime) (util/javatime-day-of-year) response) 
       (-> 0 response))))
  ([pixel_map pixel_models query-day]
   (let [values (map #(time-of-change % query-day (get pixel_map "px") (get pixel_map "py")) pixel_models)]
     (last (sort-by :val values)))))

(defn time-since-change
  "Return cumulative distance to previous break"
  ([model query-day x y]
   (let [change-prob (get model "chprob")
         break-day   (-> (get model "bday") (util/to-ordinal)) 
         query-ord   (-> query-day (util/to-ordinal))
         distance    (if (= 1.0 change-prob) (- query-ord break-day) 0)] ; can't use nil, dont think 0 is appropriate
     (hash-map :pixelx x :pixely y :val distance)))
  ([pixel_map pixel_models query-day]
   (let [values (map #(time-since-change % query-day (get pixel_map "px") (get pixel_map "py")) pixel_models)]
     (last (filter (fn [i] (some? (:val i))) (sort-by :val values))))))

(defn magnitude-of-change
  "Return severity of spectral shift"
  ([model query-day x y]
   (let [change-prob (get model "chprob")
         break-day   (-> (get model "bday") (util/to-ordinal)) 
         query-year  (-> query-day (util/to-javatime) (util/javatime-year))
         break-year  (-> break-day (util/ordinal-to-javatime) (util/javatime-year))
         magnitudes  [(get model "grmag") (get model "remag") (get model "nimag") (get model "s1mag") (get model "s2mag")]
         euc-norm    (math/sqrt (reduce + (map #(math/expt % 2) magnitudes)))
         response    #(hash-map :pixelx x :pixely y :val %)]
     (if (= true (= query-year break-year) (= 1.0 change-prob))
       (-> euc-norm (response))
       (-> 0 (response)))))
  ([pixel_map pixel_models query-day]
   (let [values (map #(magnitude-of-change % query-day (get pixel_map "px") (get pixel_map "py")) pixel_models)]
     (last (sort-by :val values)))))

(defn length-of-segment
  "Return length of change segment in days"
  ([model query-day x y]
   (let [query-ord (-> query-day (util/to-ordinal))
         start-day (-> (get model "sday") (util/to-ordinal))
         end-day   (-> (get model "eday") (util/to-ordinal))
         startends [(- query-ord start-day) (- query-ord end-day)]
         positives (filter (fn [i] (> i 0)) startends)
         minimum   (if (= 0 (count positives)) 0 (apply min positives))]
     (hash-map :pixelx x :pixely y :val minimum)))
  ([pixel_map pixel_models query-day]
   (let [values (map #(length-of-segment % query-day (get pixel_map "px") (get pixel_map "py")) pixel_models)]
     (first (sort-by :val values)))))

(defn curve-fit
  "Return Curve QA for point in time"
  ([model query-day x y]
   (let [query-ord (-> query-day (util/to-ordinal))
         curve-qa  (get model "curqa")
         start-day (-> (get model "sday") (util/to-ordinal)) 
         end-day   (-> (get model "eday") (util/to-ordinal)) 
         value     (if (<= start-day query-ord end-day) curve-qa 0)]
     (hash-map :pixelx x :pixely y :val value)))
  ([pixel_map pixel_models query-day]
   (let [values (map #(curve-fit % query-day (get pixel_map "px") (get pixel_map "py")) pixel_models)]
     (last (sort-by :val values)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;    CLASSIFICATION PRODUCTS    ;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn ismap?
  "Returns boolean true / false if input is a map "
  [input]
  (= (type input) clojure.lang.PersistentArrayMap))

(defn falls_between
  "Used to reduce a sorted list of maps to the members
  surrounding a query day"
  [mapA mapB]
  (if (ismap? mapA)
      (if (= true (:follows mapA) (:precedes mapB))
          (do [mapA mapB])
          (do mapB))       
      (do mapA)))

(defn model-class
  "Return the index of the desired classification confidence"
  [model query-day rank]
  (let [probs (get model "probs")
        sorted (reverse (sort probs)) 
        position (.indexOf probs (nth sorted rank))
        sday (-> (get model "sday") (util/to-ordinal)) 
        eday (-> (get model "eday") (util/to-ordinal)) 
        intersects (<= sday query-day eday)
        precedes   (< query-day sday)
        follows    (> eday query-day)]
    (hash-map :intersects intersects
              :precedes   precedes
              :follows    follows
              :value (nth (:lc_map config) position))))

(defn landcover
  [pixel_map models query_day rank]
  (let [px (get pixel_map "px") ; do we need these?
        py (get pixel_map "py")
        sorted_models     (sort-by-sday models)
        query_ord         (-> query_day (util/to-ordinal))
        first_start_day   (-> (first sorted_models) (get "sday") (util/to-ordinal))
        last_end_day      (-> (last sorted_models)  (get "eday") (util/to-ordinal))
        model_classes     (map #(model-class % query_ord rank) sorted_models)
        intercepted_model (first (filter :intersects model_classes))
        fell_between      (reduce falls_between model_classes)]

    ; if query_ord is less than first_start_day 
    ;
    ;    if (:fill_begin config)
    ;       return (model-class (first sorted_models) query_ord "first")
    ;    return (:lc_insuff (:lc_defaults config))


    ; if query_ord is greater than the last_end_day

    ;; (if intercepted_model
    ;;   (:value intercepted_model))

    (cond
      (= true (< query_ord first_start_day)  (:fill_begin config)) 
        (:value (first model_classes))
      (= true (< query_ord first_start_day)) 
        (:lc_insuff (:lc_defaults config))
      (= true (> last_end_day query_ord) (:fill_end config)) 
        (:value (last model_classes))
      (= true (> last_end_day query_ord))    
        (:lc_insuff (:lc_defaults config))
      (not (nil? intercepted_model))         
        (:value intercepted_model)
      (= true (:fill_samelc config) (= (:value (first fell_between)) (:value (last fell_between))))
        (:value (last fell_between))
      (= true (:fill_difflc config) )
      ) ;close cond






    )
)






;; (defn primary-landcover-orig
;;   "Return highest landcover class value for intersecting segments, between segments
;;   or outside time series"
;;   ([model query-day x y]
;;    (let [start-day (-> (get model "sday") (util/to-ordinal)) 
;;          end-day   (-> (get model "eday") (util/to-ordinal)) 
;;          query-ord (-> query-day (util/to-javatime) (util/javatime-to-ordinal))
;;          query-inclusive (<= start-day query-ord end-day)
;;          max-prob (max (get model "clprob"))
;;          between_seg_val (get model "betweensegmentvalue")]
;;      (hash-map :pixelx x :pixely y )
;;      )

;;   ) 
;;   ([pixel_map pixel_models query-day]
;;    (let [sorted_models (sort-by-sday pixel_models)
;;          values (map #(primary-landcover % query-day (get pixel_map "px") (get pixel_map "py")) sorted_models)]

;;      )))

(defn secondary-landcover
  "Return the second highest landcover class value"
  ([model query-day x y])
  ([pixel_map pixel_models query-day])
)

(defn primary-landcover-confidence
  "Return the landcover probability for the highest landcover class value"
  ([model query-day x y])
  ([pixel_map pixel_models query-day])
)

(defn secondary-landcover-confidence
  "Return the landcover probability for the 2nd highest landcover class value"
  ([model query-day x y])
  ([pixel_map pixel_models query-day])
)

(defn variable-juxt
  "Return a juxt function that returns the values for the specified map keys"
  [mapkeys]
  (apply juxt (map (fn [i] #(get % i)) mapkeys)))

(defn merge-maps-by-keys
  "Merge two lists of hash-maps, joining the list members by the key values
   specified in the mapkeys argument"
  [maplist1 maplist2 mapkeys]
  (let [conc_lists (concat maplist1 maplist2)
        grouped_lists (group-by (variable-juxt mapkeys) conc_lists)]
    (map #(merge (first %) (last %)) (vals grouped_lists))))

(defn data
  "Returns a flat list of product values from JSON of a chips worth of CCDC results"
  [segments_json predictions_json product_type queryday]
  (let [; merge segments and predictions by px, py, cx, cy, sday and eday
        segments_predictions (merge-maps-by-keys segments_json predictions_json ["px" "py" "cx" "cy" "sday" "eday"])
        ; group segments by pixel coordinates
        pixel_segments (util/coll-groups segments_predictions ["px" "py"])
        ; map the products function across the pixel segments. Returns a flat
        ; collection, one hash map per pixel coordinate pair.
        product_fn (-> (str "lcmap.gaia.products/" product_type) (symbol) (resolve))
        pixel_array (map #(product_fn (first %) (last %) queryday) pixel_segments)
        ; group product coll by row 
        ; [{:pixely 3159045} [{:pixely 3159045, :pixelx -2114775, :val 6290},...] ...]
        row_groups (util/coll-groups pixel_array [:pixely]) 
        ; sort row group values by pixelx ascending 
        sort-pixelx-fn (fn [i] (hash-map (:pixely (first i)) (sort-by :pixelx (last i))))
        sorted-x-vals (map sort-pixelx-fn row_groups)
        ; sort the rows by the pixely key ascending
        sorted-y-rows (sort-by (fn [i] (first (keys i))) sorted-x-vals)]
    ; finally, flatten to a one dimensional list
    (util/flatten-vals sorted-y-rows :val)))

