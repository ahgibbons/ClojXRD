;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Math Functions
;;;
;;;;;;;;;;;;;;;;;;;;;

(def pi Math/PI)

(defn rad->deg [x]
  (* x (/ 180 pi))
)

(defn deg->rad [x]
  (* x (/ pi 180)))

(defn sq [x] (* x x))

(defn factorial [n]
  (loop [n' n
         acc 1]
    (if (zero? n')
      acc
      (recur (dec n') (* acc n'))))
)

(defn exp [x n]
  (loop [acc 1 n' n]
    (if (zero? n') acc
        (recur (* x acc) (dec n'))))
)

(defn v-add [& vectors]
  (reduce (fn [v1 v2] (map + v1 v2)) vectors))

(defn vec-mult [v a]
  (map #(* a %) v))

(defn vec-div [v a]
  (map #(/ % a) v))


(defn v-cross [v1 v2]
  (let [x1 (nth v1 0)
        y1 (nth v1 1)
        z1 (nth v1 2)
        x2 (nth v2 0)
        y2 (nth v2 1)
        z2 (nth v2 2)]
    [(- (* y1 z2) (* y2 z1))
     (- (* x2 z1) (* x1 z2))
     (- (* x1 y2) (* x2 y1))]))

(defn v-dot [v1 v2]
  (reduce + (map * v1 v2)))

(defn v-mag-sq [x]
  (reduce + (map * x x)))

(defn v-mag [x]
  (Math/sqrt (v-mag-sq x)))

(defn volume [a1 a2 a3]
  (v-dot a1 (v-cross a2 a3)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Crystal Stuff
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;




;; Lamda Cu
(def lamda-cu 1.54059)

(defn vectors->lattice-params [a1 a2 a3]
  (let [a (v-mag a1)
        b (v-mag a2)
        c (v-mag a3)
        alpha (rad->deg (Math/acos (/ (v-dot a2 a3) (* b c))))
        beta (rad->deg (Math/acos (/ (v-dot a1 a3) (* a c))))
        gamma (rad->deg (Math/acos (/ (v-dot a1 a2) (* b c))))]
    {:a a :b b :c c :alpha alpha :beta beta :gamma gamma})
)


(defn reciprocal-vector [a1 a2 a3]
  (vec-div (v-cross a2 a3) (volume a1 a2 a3)))

(defn g-hkl [a1 a2 a3 h k l]
  (let [b1 (reciprocal-vector a1 a2 a3)
        b2 (reciprocal-vector a2 a3 a1)
        b3 (reciprocal-vector a3 a1 a2)]
    (v-add (vec-mult b1 h) (vec-mult b2 k) (vec-mult b3 l))))

(defn d-spacing [a1 a2 a3 h k l]
  (/ 1 
     (v-mag (g-hkl a1 a2 a3 h k l))))

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




(defn d-spacing-cubic [a h k l] 
  (/ a (Math/sqrt (+ (* h h) (* k k) (* l l))))
)


(defn d-spacing-ortho [a b c h k l]
  (/ 1 (Math/sqrt 
        (+ (sq (/ h a)) (sq (/ k b)) (sq (/ l c))))))


(defn multiplicity-help [nzeros nequal]
  (/ 
      (* (factorial 3) (exp 2 (- 3 nzeros)))
      (factorial nequal))
)

(defn multiplicity [h k l]
 (let [nzeros (count (filter zero? (list h k l)))
       nequal (case (count (distinct (list h k l)))
                1 3
                2 2
                3 1)]
   (multiplicity-help nzeros nequal))
)

(def hkls '((1 0 0) (1 1 0) (1 1 1) (2 0 0) (2 1 0) (2 2 0) (2 2 2) (2 1 1) (2 2 1)))

(defn bragg-angle
  [lambda d ]
  (rad->deg (Math/asin (/ lambda (* 2 d))))
)

(defn cartesian-prod [colls]
  (if (empty? colls)
    '(())
    (for [more (cartesian-prod (rest colls))
          x (first colls)]
      (cons x more))))

(defn .+ [x1 x2]
  {:real (+ (:real x1) (:real x2)) :imag (+ (:imag x1) (:imag x2))})


(defn structure-factor [f h k l x y z]
  (let [x (*  -2 pi (+ (* h x) (* k y) (* l z)))
       real (* f (Math/cos x))
       imag (* f (Math/sin x))]
  {:real real :imag imag}))

(def Na-scat-const {:a1 4.7626 :b1 3.285 :a2 3.1736 :b2 8.8422
                    :a3 1.2674 :b3 0.3136 :a4 1.1128 :b4 129.424 :c 0.676})

;; Na 4.7626 3.285  3.1736  8.8422  1.2674  0.3136  1.1128  129.424 0.676 
;; Na1+ 3.2565 2.6671 3.9362 6.1153 1.3998  0.2001  1.0032  14.039  0.404 


(defn atom-scattering-factor [atom-scat-consts G]
  (let [a1 (atom-scat-consts :a1)
        b1 (:b1 atom-scat-consts)
        a2 (:a2 atom-scat-consts)
        b2 (:b2 atom-scat-consts)
        a3 (:a3 atom-scat-consts)
        b3 (:b3 atom-scat-consts)
        a4 (:a4 atom-scat-consts)
        b4 (:b4 atom-scat-consts)
        c (:c atom-scat-consts)
        fterm (fn [a b] (* a (Math/exp (/ (* -1 b G) (* 4 pi)))))]
    (+ (fterm a1 b1) (fterm a2 b2) (fterm a3 b3) (fterm a4 b4) c)
  ))

(defn structure-factor-atom [ato h k l]
  (let [f (:f ato)
        x (:x ato)
        y (:y ato)
        z (:z ato)]
    (structure-factor f h k l x y z)))


(defn F-hkl [atoms h k l]
   (reduce .+ (map (fn [ato] (structure-factor-atom ato h k l)) atoms)
  ))

(def atom1 {:f 1.0 :x 0 :y 0 :z 0})
(def atom2 {:f 1.0 :x 0.5 :y 0.5 :z 0.5})

;(defn x-ray-reflection 
;  [lambda a hkl]
;  (let [h (nth hkl 0)
;        k (nth hkl 1)
;        l (nth hkl 2)
;        d (d-spacing-cubic a h k l)
;        theta (bragg-angle lambda d)]
;    (list (* 2 theta) (multiplicity h k l)))
;)

(def t-a1 [5 0 0])
(def t-a2 [0 5 0])
(def t-a3 [0 0 5])


(def Crystal-System {:P23 :cubic
                     :F23 :cubic
                     :I23 :cubic
                     :P213 :cubic
                     :I213 :cubic
                     :Pm-3m :cubic
                     :Pn-3 :cubic
                     :Fm-3 :cubic
                     :Fd-3 :cubic
                     :Im-3 :cubic
                     :Pa-3 :cubic
                     :Ia-3 :cubic
                     :P432 :cubic
                     :P4232 :cubic
                     :F432 :cubic
                     })

(defn lattice-vectors-from-params [unit-cell]
  (case (Crystal-System (unit-cell :space-group)) 
    :cubic (let [a (unit-cell :a)]
             {:a [a 0 0]
              :b [0 a 0]
              :c [0 0 a]})))


(def simple-Na {:a 3.5 :b 3.5 :c 3.5 :alpha 90 :beta 90 :gamma 90 :space-group :Pm-3m
                :pos [{:atom :Na :x 0.0 :y 0.0 :z 0.0}]})
(def simple-a [3.5 0 0])
(def simple-b [0 3.5 0])
(def simple-c [0 0 3.5])

(def simple-G  (g-hkl simple-a simple-b simple-c 1 0 0 ))

(def NaCl {:a 5.58813 :b 5.58813 :c 5.58813 :alpha 90 :beta 90 :gamma 90 :space-group :P1 
            :pos [{:atom :Na1+ :x 0.0 :y 0.0 :z 0.0}
                  {:atom :Na1+ :x 0.0 :y 0.5 :z 0.5}
                  {:atom :Na1+ :x 0.5 :y 0.0 :z 0.5}
                  {:atom :Na1+ :x 0.5 :y 0.5 :z 0.0}
                  {:atom :Cl1- :x 0.0 :y 0.0 :z 0.5}
                  {:atom :Cl1- :x 0.0 :y 0.5 :z 0.0}
                  {:atom :Cl1- :x 0.5 :y 0.0 :z 0.0}
                  {:atom :Cl1- :x 0.5 :y 0.5 :z 0.5}]})

(defn XRD [lambda-xray unit-cell ])
