
(def pi Math/PI)
(def lambda0 0.154) ;; nanometers
(defn deg->rad [x] (* x (/ pi 180.0)))
(defn rad->deg [x] (* x (/ 180 pi)))
(def min-angle 1)


(defn sind [x]
    (Math/sin (deg->rad x)))

(defn cosd [x]
    (Math/cos (deg->rad x)))

(defn acosd [x]
    (rad->deg (Math/acos x)))

(defn asind [x]
    (rad->deg (Math/asin x)))
;; Vector Operations

(defn v* [a xs] (map #(* % a) xs))
(defn dot [v1 v2] (reduce + (map * v1 v2)))
(defn cross [v1 v2] (
        let [
             x1 (nth v1 0) x2 (nth v2 0)
             y1 (nth v1 1) y2 (nth v2 1)
             z1 (nth v1 2) z2 (nth v2 2)]
        [(- (* y1 z2) (* y2 z1))
          (- (* z1 x2) (* x1 z2))
          (- (* x1 y2) (* y1 x2))]))



(defn vector-volume [a b c]
    (dot a (cross b c)))

;;; Unit Cell Calculations

(defn cell-volume [a b c alpha beta gamma] 
    (let [alpha-r (deg->rad alpha)
          beta-r  (deg->rad beta)
          gamma-r (deg->rad gamma)
          cos-a   (Math/cos alpha-r)
          cos-b   (Math/cos beta-r)
          cos-g   (Math/cos gamma-r)
          ]
        (* a b c 
           (Math/sqrt 
               (+ 1 
                  (* 2 (* cos-a cos-b cos-g)) 
                  (- (+ (* cos-a cos-a) (* cos-b cos-b) (* cos-g cos-g))))))))


(defn d-spacing [a b c alpha beta gamma h k l]
    (let [V (cell-volume a b c alpha beta gamma)
          alpha-r (deg->rad alpha)
          beta-r (deg->rad beta)
          gamma-r (deg->rad gamma)
          sin-a (Math/sin alpha-r)
          sin-b (Math/sin beta-r)
          sin-g (Math/sin gamma-r)
          cos-a (Math/cos alpha-r)
          cos-b (Math/cos beta-r)
          cos-g (Math/cos gamma-r)
          t1 (* h h b b c c sin-a sin-a)
          t2 (* k k a a c c sin-b sin-b)
          t3 (* l l a a b b sin-g sin-g)
          t4 (* h k a b c c (* (- (* cos-a cos-b) cos-g)))
          t5 (* h l a b b c (* (- (* cos-a cos-g) cos-b)))
          t6 (* l k a a b c (* (- (* cos-a cos-g) cos-a)))

          td (+ t1 t2 t3 t4 t5 t6)
          tn (* V V)]
        (Math/sqrt (/ tn td))
    )
)

(defn dspacing
    [h k l a]
    (/ a (Math/sqrt (+ (* h h) (* k k) (* l l)))))

(defn dspacing2
    [hkl a]
    (dspacing (nth hkl 0) (nth hkl 1) (nth hkl 2) a))

(defn min-spacing
    [wavelength twotheta]
    (/ wavelength (* 2 (Math/sin  (deg->rad (/ twotheta 2))))))

(defn valid-hkl? [hkl a min-d]
    (and (> (dspacing2 hkl a) min-d)
         (not (= (nth hkl 0) (nth hkl 1) (nth hkl 2)))))

(defn valid-hkl?-alt [h k l a min-d]
    (and 
         (not (= h k l 0))
         (> (dspacing h k l a) min-d))    
)

(defn hkl-values
    [a min-d]
    (let [hmax (Math/floor (/ a min-d)) 
          kmax (Math/floor (/ a min-d)) 
          lmax (Math/floor (/ a min-d))
          hkls (for [h (range (- hmax) hmax) k (range (- kmax) kmax) l (range (- lmax) lmax)]
                        [(int h) (int k) (int l)])]
            (filter (fn [h k l] (valid-hkl?-alt h k l a min-d)) hkls
            )))


;;   
;;  1. define unit cell
;;  2. generate miller indices
;;  3. Calculate d-spacing
;;  4. Apply Braggs law 2dsina = lambda
;;  5. calculate structure factor F
;;

;;;; Unit Cell

(def NaCl [{:a 5.58813 :b 5.58813 :c 5.58813 :alpha 90 :beta 90 :gamma 90 :space-group :P1 
            :pos [{:atom :Na1+ :x 0.0 :y 0.0 :z 0.0}
                  {:atom :Na1+ :x 0.0 :y 0.5 :z 0.5}
                  {:atom :Na1+ :x 0.5 :y 0.0 :z 0.5}
                  {:atom :Na1+ :x 0.5 :y 0.5 :z 0.0}
                  {:atom :Cl1- :x 0.0 :y 0.0 :z 0.5}
                  {:atom :Cl1- :x 0.0 :y 0.5 :z 0.0}
                  {:atom :Cl1- :x 0.5 :y 0.0 :z 0.0}
                  {:atom :Cl1- :x 0.5 :y 0.5 :z 0.5}]}])

(def NaCl-hkl (hkl-values (get NaCl :a) (min-spacing lambda0 min-angle)))


;(defn reciprocal-lattice-vec 
;    [a b c alpha beta gamma]
;    (let [fac (/ (* 2 pi ) (determinant a b c))]
;        [ (v* fac (cross b c )) (v* fac (cross c a)) (v* fac (cross a b))])
;    )



(def atoms {:Na1+   {  :a1	 3.2565 	:b1 2.6671 	:a2 3.9362 	:b2 6.1153 	:a3 1.3998 	:b3 0.2001 	:a4 1.0032 	:b4 14.039 	:c 0.404}
            :Cl1-   { :a1 18.2915 :b1 0.0066 :a2 7.2084 :b2 1.1717 :a3 6.5337 :b3 19.5424 :a4 2.3386 :b4 60.4486 :c -16.378}
            
            })   
;`(Cl1- 	 18.2915 	 0.0066 	 7.2084 	 1.1717 	 6.5337 	 19.5424 	 2.3386 	 60.4486 	 -16.378) 
