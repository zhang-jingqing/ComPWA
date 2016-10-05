;;;*************************
;;;* DECAY VIOLATION RULES *
;;;*************************

(defrule remove-bad-L
	(declare (salience 99))
	?decay <- (Decay 
				(quantum_number_name ?qn_name) 
				(mother ?mother_id) 
				(daughters ?daughter1_id ?daughter2_id $?others)
				(required_variable_names "angular-momentum" $?other_rvns) 
				(violating_quantum_number_list $?violating_quantum_number_list)
			  )
	(test (not (= 0 (str-compare ?qn_name "spinwave"))))
	=>
	(bind ?angular_momentum (get-spin-qn-with-unique-id (get-required-variable "angular-momentum" ?decay)))
	; in the helicity formalism we throw out all z components of L except 0
	(if (<> 0 (fact-slot-value ?angular_momentum z_component_numerator))
	then
		;(printout t "removing " ?decay crlf)
		(retract ?decay)
    )
)

(defrule check-spin
	(declare (salience 99))
	(SpinQuantumNumber (unique_id ?mother_id) (numerator ?num_mother) (denominator ?denom_mother))
	(SpinQuantumNumber (unique_id ?daughter1_id) (numerator ?num_daughter1) (denominator ?denom_daughter1))
	(SpinQuantumNumber (unique_id ?daughter2_id) (numerator ?num_daughter2) (denominator ?denom_daughter2))
	?decay <- (Decay 
				(quantum_number_name "spin")
				(mother ?mother_id)
				(daughters ?daughter1_id ?daughter2_id $?others)
				(required_variable_names "angular-momentum" $?other_rvns)
				(violating_quantum_number_list $?violating_quantum_number_list)
			  )
	(test (not (member$ "spin" ?violating_quantum_number_list)))
	=>
	;get the required information
	(bind ?angular_momentum (get-spin-qn-with-unique-id (get-required-variable "angular-momentum" ?decay)))
	(bind ?L (/ (fact-slot-value ?angular_momentum numerator) (fact-slot-value ?angular_momentum denominator)))
	
	
	(bind ?violated TRUE)
	(bind ?comb1 (create$))
	(loop-for-count
		(?i
			(abs (- (/ ?num_daughter1 ?denom_daughter1) (/ ?num_daughter2 ?denom_daughter2)))
			(+ (/ ?num_daughter1 ?denom_daughter1) (/ ?num_daughter2 ?denom_daughter2))
		)
	do
		(bind ?comb1 (insert$ ?comb1 1 ?i))
	)
	;(printout t ?L " " ?comb1 crlf)
	;(printout t ?L " " (/ ?num_mother ?denom_mother) " " (/ ?num_daughter1 ?denom_daughter1) " "  (/ ?num_daughter2 ?denom_daughter2) crlf)
	(foreach ?val ?comb1
	do
		(if (and
				(<= (/ ?num_mother ?denom_mother) 
					(+ ?val ?L)
				)
				(>= (/ ?num_mother ?denom_mother) 
					(abs (- ?val ?L))
				)
			)
		then
			(bind ?violated FALSE)
			(break)
		)
	)
	
	(if ?violated
	then
		(if (is-qn-conserved "spin")
		then
			;(printout t "decay violates angular momentum conservation!" crlf)
			(retract ?decay)
	  	else
	  		(modify ?decay (violating_quantum_number_list ?violating_quantum_number_list "spin"))
	  	)
	)
)

(defrule check-helicity
	(declare (salience 99))
	(SpinQuantumNumber (unique_id ?mother_id) (numerator ?num_mother) (denominator ?denom_mother))
	(SpinQuantumNumber (unique_id ?daughter1_id) (z_component_numerator ?z_num_daughter1) (denominator ?denom_daughter1))
	(SpinQuantumNumber (unique_id ?daughter2_id) (z_component_numerator ?z_num_daughter2) (denominator ?denom_daughter2))
	?decay <- (Decay (quantum_number_name "spin") (mother ?mother_id) (daughters ?daughter1_id ?daughter2_id $?others))
	=>
	(if (< (/ ?num_mother ?denom_mother) (abs (- (/ ?z_num_daughter1 ?denom_daughter1) (/ ?z_num_daughter2 ?denom_daughter2))))
	then
	  ;(printout t "decay violates angular momentum conservation!" crlf)
	  (retract ?decay)
	)
)

(defrule check-parity-helicity
	(declare (salience 99))

	?decay_tree <- (DecayTree (decays $?decays) (available_waves $?available_waves))
	
	=>
	;(printout t "checking for parity violating decays in helicity basis (only A00 amplitudes)" crlf)

    (if (is-qn-conserved "parity")
    then
	
		(foreach ?decay ?decays
			(bind ?mother_index (fact-slot-value ?decay mother))
			(bind ?d1_index (nth$ 1 (fact-slot-value ?decay daughters)))
			(bind ?d2_index (nth$ 2 (fact-slot-value ?decay daughters)))
				
			; get all info that is needed
			(bind ?parity_mother (get-qn-value ?mother_index ?decay_tree "parity"))
			(bind ?parity_daughter1 (get-qn-value ?d1_index ?decay_tree "parity"))
			(bind ?parity_daughter2 (get-qn-value ?d2_index ?decay_tree "parity"))
	
			(bind ?spin_mother (get-spin-qn-with-unique-id (get-qn-value ?mother_index ?decay_tree "spin")))
			(bind ?spin_daughter1 (get-spin-qn-with-unique-id (get-qn-value ?d1_index ?decay_tree "spin")))
			(bind ?spin_daughter2 (get-spin-qn-with-unique-id (get-qn-value ?d2_index ?decay_tree "spin")))
	
			(if (and
					(= (fact-slot-value ?spin_daughter1 z_component_numerator) 
						0
					)
					(= (fact-slot-value ?spin_daughter2 z_component_numerator) 
						0
					)
				)
			then
	    		(bind ?intrinsic_parity_product (* ?parity_mother (* ?parity_daughter1 ?parity_daughter2)))
	  		    (bind ?spin_difference 
	    			(abs
	    				(- 
	    					(+ 
	    						(/ (fact-slot-value ?spin_daughter1 numerator) (fact-slot-value ?spin_daughter1 denominator)) 
	    						(/ (fact-slot-value ?spin_daughter2 numerator) (fact-slot-value ?spin_daughter2 denominator))
	    					)
	    					(/ (fact-slot-value ?spin_mother numerator) (fact-slot-value ?spin_mother denominator))
	    				)
	    			)
	    		)
	    		(bind ?spin_factor_part -1)
	    		(if (= 0 (mod ?spin_difference 2)) then (bind ?spin_factor_part 1))
				
				;(printout t "checking if " ?intrinsic_parity_product " " ?spin_factor_part crlf)
				(if (= -1 (* ?intrinsic_parity_product ?spin_factor_part))
				then
					;(printout t "decay violates parity!" crlf)
					(if (is-qn-conserved "parity")
					then
						(retract ?decay_tree)
	  					(break)
	  				else
	  					(bind ?violating_quantum_number_list (fact-slot-value ?decay violating_quantum_number_list))
	  					(if (not (member$ "parity" ?violating_quantum_number_list))
	  					then
	  						(modify ?decay (violating_quantum_number_list ?violating_quantum_number_list "parity"))
	  					)
	  				)
				)
			)
		)
	)
)
