;2D PC. Triangular air holes drilled in a dielectric substrate
;Constants
;(define r 0.48)
(define r 0.4)
(define a1 (vector3 1 0 0))
(define a2 (vector3 (/ 1 2) (/ (sqrt 3) 2) 0))
(define a3 (vector3 0 0 1))

;Geometry
;Default material to fill space.
(set! default-material (make dielectric (epsilon 13)))
;Lattice
(set! geometry-lattice (make lattice (basis1 a1) (basis2 a2) (basis3 a3) (basis-size 1 1 1) (size 1 1 no-size) ))

;Objects
(set! geometry (list
		(make cylinder (material air) (center 0 0 0) (radius r) (height no-size) (axis 0 0 1))
))

(set! mesh-size 20)
;Parameters
;Irreducible B.Z. for triangular lattice
(define Gamma (vector3 0 0 0))
(define M (vector3 0 0.5 0))
(define K (vector3 (/ 1 -3) (/ 1 3) 0))
;(define M (vector3 (/ 2 3) (/ 1 3) 0))
(define interp 20)
(define res 64)

(set! k-points (interpolate interp(list Gamma M K Gamma)))
(set! resolution res)
(set! num-bands 4)
(set! filename-prefix false)
(run-te fix-efield-phase output-efield display-group-velocities)
(run-tm fix-efield-phase output-efield display-group-velocities)
