;;1D PC. GaAs Bulk. 
;;Geometry File
;;Constants
(define GaAs (make dielectric (epsilon 13)))

;;PUC geometry
(set! geometry (list(
		     make block (material GaAs) (center 0 0 0) (size 1 no-size no-size) (e1 1 0 0) (e2 0 1 0) (e3 0 0 1) )
))

(set! geometry-lattice (
			make lattice (basis1 1 0 0) (basis2 0 1 0) (basis3 0 0 1) (basis-size 1 1 1) (size 1 no-size no-size)
))

(set! resolution (vector3 32 1 1))
(set! mesh-size 3)

;GaAs Bulk

;Parameters.
(define Gamma (vector3 0 0 0))
(define X (vector3 0.5 0 0))

;Irreducible B.Z.
(set! k-points (interpolate 20 (list Gamma X)))

(set! num-bands 4)
(set! filename-prefix false)
(run fix-efield-phase output-efield)
