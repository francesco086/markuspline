
MODULE example
	IMPLICIT NONE
	

	CONTAINS

	FUNCTION trial_function(i,x)
		IMPLICIT NONE
		REAL(KIND=8) :: trial_function
		INTEGER, INTENT(IN) :: i
		REAL(KIND=8), INTENT(IN) :: x
		
		!SELECT CASE(i)
		!CASE(0)
		!	trial_function=DSQRT(1.d0+x)
		!CASE(1)
		!	trial_function=1.d0/(2.d0*((1.d0+x)**(1./2.)))
		!CASE(2)
		!	trial_function=-1.d0/(4.d0*((1.d0+x)**(3./2.)))
		!CASE(3)
		!	trial_function=3.d0/(8.d0*((1.d0+x)**(5./2.)))
		!CASE(4)
		!	trial_function=-15.d0/(16.d0*((1.d0+x)**(7./2.)))
		!CASE(5)
		!	trial_function=105.d0/(32.d0*((1.d0+x)**(9./2.)))
		!END SELECT

		SELECT CASE(i)
		CASE(0)
			trial_function=1.d0/(1.d0+x)
		CASE(1)
			trial_function=-1.d0/((1.d0+x)**2)
		CASE(2)
			trial_function=2.d0/((1.d0+x)**3)
		CASE(3)
			trial_function=-6.d0/((1.d0+x)**4)
		CASE(4)
			trial_function=24.d0/((1.d0+x)**5)
		CASE(5)
			trial_function=-120.d0/((1.d0+x)**6)
		END SELECT

	END FUNCTION trial_function


END MODULE example
