;GaAs_Air 1D PC.

(define GaAs (make dielectric (epsilon 13)))
;(define GaAs2 (make dielectric (epsilon 10)))

(set! geometry-lattice (
			make lattice (basis-size 1 1 1) (basis1 1 0 0) (basis2 0 1 0) (basis3 0 0 1) (size 1 no-size no-size)
))

(set! geometry (list
		      ; MPB creates block's span as (center-size/2, center+size/2). So the following two lines creates the blocks
	              ; GaAs:(-0.25, 0.25) Air:(0.25, 0.75) (-0.25,0) is translated to (0.75,1)
		      ; (make block (material GaAs) (center 0 0 0) (size 0.5 no-size no-size) (e1 1 0 0) (e2 0 1 0) (e3 0 0 1))
		      ; (make block (material air) (center 0.5 0 0) (size 0.5 no-size no-size) (e1 1 0 0) (e2 0 1 0) (e3 0 0 1))

		      ; (make block (material GaAs) (center 1.0 0 0) (size 0.5 no-size no-size) (e1 1 0 0) (e2 0 1 0) (e3 0 0 1))
	              ; (make block (material air) (center 1.5 0 0) (size 0.5 no-size no-size) (e1 1 0 0) (e2 0 1 0) (e3 0 0 1))

                      (make block (material GaAs) (center 0.25 0 0) (size 0.5 no-size no-size))
		      (make block (material air) (center 0.75 0 0) (size 0.5 no-size no-size))
))

(set! resolution (vector3 64 64 64)) ;resolution in each basis direction.
(set! mesh-size 10)
;Parameters set for computation.

(define Gamma (vector3 0 0 0))
(define X (vector3 0.5 0 0))

;Irreducible B.Z
(set! k-points (interpolate 20 (list Gamma X)))

;Number of bands
(set! num-bands 4)
;(set! num-bands 2)
(set! filename-prefix false)
(run fix-efield-phase output-efield display-group-velocities)
