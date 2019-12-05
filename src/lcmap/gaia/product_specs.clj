(ns lcmap.gaia.product-specs
  (:require [clojure.spec.alpha :as spec]
            [clojure.walk :refer [keywordize-keys]]))

(spec/def ::intr float?)
(spec/def ::magr float?)
(spec/def ::rmsr float?)

(spec/def ::count_seven #(= 7 (count %)))
(spec/def ::doubles (spec/coll-of double?))
(spec/def ::coefr (spec/and ::doubles ::count_seven))

; too restrictive, and dumb
;(spec/def ::coord_bounds #(or (< -9999999 % -1000000) (> 9999999 % 1000000)))
;(spec/def ::coord (spec/and integer? ::coord_bounds))
(spec/def ::coord integer?)

(spec/def ::not_default #(not (re-matches #"0001-01-01" %)))
(spec/def ::day_fmt #(some? (re-matches #"[0-9]{4}-[0-9]{2}-[0-9]{2}" %)))
(spec/def ::days (spec/and string? ::day_fmt ::not_default))

; segments predicates
(spec/def ::px ::coord)
(spec/def ::py ::coord)
(spec/def ::sday ::days)
(spec/def ::eday ::days)
(spec/def ::bday ::days)
(spec/def ::chprob float?)
(spec/def ::curqa integer?)
(spec/def ::blint ::intr)
(spec/def ::grint ::intr)
(spec/def ::niint ::intr)
(spec/def ::reint ::intr)
(spec/def ::s1int ::intr)
(spec/def ::s2int ::intr)
(spec/def ::thint ::intr)
(spec/def ::blmag ::magr)
(spec/def ::grmag ::magr)
(spec/def ::nimag ::magr)
(spec/def ::remag ::magr)
(spec/def ::s1mag ::magr)
(spec/def ::s2mag ::magr)
(spec/def ::thmag ::magr)
(spec/def ::blrmse ::rmsr)
(spec/def ::grrmse ::rmsr)
(spec/def ::nirmse ::rmsr)
(spec/def ::rermse ::rmsr)
(spec/def ::s1rmse ::rmsr)
(spec/def ::s2rmse ::rmsr)
(spec/def ::thrmse ::rmsr)
(spec/def ::blcoef ::coefr)
(spec/def ::grcoef ::coefr)
(spec/def ::nicoef ::coefr)
(spec/def ::recoef ::coefr)
(spec/def ::s1coef ::coefr)
(spec/def ::s2coef ::coefr)
(spec/def ::thcoef ::coefr)
(spec/def ::segment (spec/keys :req-un [::px ::py ::sday ::eday ::bday ::chprob ::curqa 
                                        ::blint  ::grint  ::niint  ::reint  ::s1int  ::s2int  ::thint
                                        ::blmag  ::grmag  ::nimag  ::remag  ::s1mag  ::s2mag  ::thmag
                                        ::blrmse ::grrmse ::nirmse ::rermse ::s1rmse ::s2rmse ::thrmse
                                        ::blcoef ::grcoef ::nicoef ::recoef ::s1coef ::s2coef ::thcoef]))
(spec/def ::segments (spec/coll-of ::segment))

; predictions predicates
(spec/def ::cx ::coord)
(spec/def ::cy ::coord)

(spec/def ::count_nine #(= 9 (count %)))
(spec/def ::prob (spec/and ::count_nine ::doubles))
(spec/def ::pday ::days)
(spec/def ::prediction (spec/keys :req-un [::cx ::cy ::px ::py ::sday ::eday ::pday ::prob]))
(spec/def ::predictions (spec/coll-of ::prediction))

(spec/def ::product_type #(contains? #{"annual-change" "curve-fit" "length-of-segment"  "magnitude-of-change" 
                                       "primary-landcover" "primary-landcover-confidence" "secondary-landcover"  
                                       "secondary-landcover-confidence" "time-of-change" "time-since-change"} %))

(spec/def ::count_ten_thousand #(= 10000 (count %)))

(defn get-spec-problems
  [explain_data]
  (first (:clojure.spec.alpha/problems explain_data)))

(defn check!
  [spec params excinfo]
  (or (some->> (spec/explain-data spec params)
               (get-spec-problems)
               (hash-map :details)
               (conj excinfo)
               (ex-info "Validation Error!")
               (throw))
      (spec/conform spec params)))

(defn segment-valid?
  [segment]
  (nil? (spec/explain-data ::segment (keywordize-keys segment))))

(defn segments-valid?
  [segments]
  (nil? (spec/explain-data ::segments (keywordize-keys segments))))

(defn prediction-valid?
  [prediction]
  (nil? (spec/explain-data ::prediction (keywordize-keys prediction))))

(defn predictions-valid?
  [predictions]
  (nil? (spec/explain-data ::predictions (keywordize-keys predictions))))

(defn segment_check
  [segment]
  (check! ::segment segment           {:message "segment-exception"      :type "data-request-error"}))

(defn prediction_check
  [prediction]
  (check! ::prediction prediction     {:message "prediction-exception"   :type "data-request-error"}))

(defn segment_coll_check
  [segments]
  (check! ::segments segments         {:message "segments-exception"     :type "data-request-error"}))

(defn prediction_coll_check
  [predictions]
  (check! ::predictions predictions   {:message "predictions-exception"  :type "data-request-error"}))

(defn product_type_check
  [product_type]
  (check! ::product_type product_type {:message "product-type-exception" :type "data-request-error"}))

(defn date_fmt_check
  [date_str]
  (check! ::days date_str             {:message "date-format-exception"  :type "data-request-error"}))

(defn output_check
  [output_values]
  (check! ::count_ten_thousand output_values {:message "output-size-exception" :type "data-generation-error"}))
