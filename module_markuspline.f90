MODULE markuspline
	IMPLICIT NONE

	LOGICAL, PRIVATE, SAVE :: MSPL_DEBUG_MODE=.TRUE.

	TYPE MSPLINE
		LOGICAL :: flag_init=.FALSE.     !Flag which specify whether the spline object has been allocated with new_MSPLINE or not
		INTEGER :: Nknots   !Number of knots minus one (the index zero is reserved for r=0)
		INTEGER :: m	     !It specifies the polynomial order, which will be 2*m+1
		INTEGER, PRIVATE :: po			!Polynomial order
		REAL(KIND=8) :: La, Lb			!They specify the range of x over which the spline has been built
		REAL(KIND=8) :: delta                      !Distance between knots,
      REAL(KIND=8), PRIVATE :: Idelta            !and its inverse (Idelta=1.d0/delta)
		REAL(KIND=8), ALLOCATABLE :: x(:)       !Grid points
		REAL(KIND=8), ALLOCATABLE :: t(:,:)       !Spline coefficients:   first index=order of derivative,   second index=nknot
		REAL(KIND=8), PRIVATE, ALLOCATABLE :: S(:,:)       !Matrix necessary to compute the spline
		REAL(KIND=8), PRIVATE, ALLOCATABLE :: ST(:,:)        !S matrix transposed, useful for improved performance
		REAL(KIND=8), PRIVATE, ALLOCATABLE :: DF(:,:)    !Factors necessary for the computation of the derivatives
      LOGICAL :: cutoff                                !Specifies whether the spline goes smoothly to zero at Lb or not
      INTEGER, PRIVATE :: fade                         !If 1 the spline vanish to zero at Lb, if 0 it does not
	END TYPE MSPLINE

   !subroutines used internally, which should not be accessible from outside
	PRIVATE :: fact
	PRIVATE :: generate_matrix_S
	PRIVATE :: generate_derivative_factors

CONTAINS

	
	FUNCTION new_MSPL(m,Nknots,La,Lb,cutoff)
	!create new spline object
		IMPLICIT NONE
		TYPE(MSPLINE) :: new_MSPL
		INTEGER, INTENT(IN) :: Nknots, m
		REAL(KIND=8), INTENT(IN) :: La, Lb
      LOGICAL, INTENT(IN) :: cutoff
		INTEGER :: i, j
		
		IF (MSPL_DEBUG_MODE) THEN
			IF (Nknots < 1) THEN
				PRINT *, "### MSPL_ERROR ### invoked by new_MSPL"
				PRINT *, "Nknots must be greater than zero"
				STOP
			END IF
			IF (m < 0) THEN
				PRINT *, "### MSPL_ERROR ### invoked by new_MSPL"
				PRINT *, "m must be equal or greater than zero"
				STOP
			END IF
         IF (La>=Lb) THEN
            PRINT *, "### MSPL_ERROR ### invoked by new_MSPL"
            PRINT *, "LB must be greater than LA"
            STOP
         END IF
		END IF
		
		new_MSPL%Nknots=Nknots
		
		new_MSPL%m=m
		new_MSPL%po=2*m+1

		new_MSPL%La=La
		new_MSPL%Lb=Lb
		new_MSPL%delta=(Lb-La)/REAL(Nknots,8)
		new_MSPL%Idelta=1.d0/new_MSPL%delta

		ALLOCATE(new_MSPL%x(-1:Nknots+1))
		DO i = -1, new_MSPL%Nknots+1, 1
			new_MSPL%x(i)=La+new_MSPL%delta*REAL(i,8)
		END DO

		ALLOCATE(new_MSPL%t(0:m,0:Nknots))
		new_MSPL%t=0.d0
		
		ALLOCATE(new_MSPL%S(0:m,0:new_MSPL%po))
		CALL generate_matrix_S(new_MSPL)
		ALLOCATE(new_MSPL%ST(0:new_MSPL%po,0:m))
		DO j = 0, m, 1
			DO i = 0, new_MSPL%po, 1
				new_MSPL%ST(i,j)=new_MSPL%S(j,i)
			END DO
		END DO

		ALLOCATE(new_MSPL%DF(0:new_MSPL%po,0:new_MSPL%po))
		CALL generate_derivative_factors(new_MSPL)

      new_MSPL%cutoff=cutoff
      IF (cutoff) THEN
         new_MSPL%fade=1
      ELSE
         new_MSPL%fade=0
      END IF
		
		new_MSPL%flag_init=.TRUE.
	
	END FUNCTION new_MSPL


	SUBROUTINE MSPL_new(m,Nknots,La,Lb,spl,cutoff)
	!create new spline object
		IMPLICIT NONE
		TYPE(MSPLINE) :: spl
		INTEGER, INTENT(IN) :: Nknots, m
		REAL(KIND=8), INTENT(IN) :: La, Lb
      LOGICAL, INTENT(IN) :: cutoff
		INTEGER :: i, j
		
		IF (MSPL_DEBUG_MODE) THEN
			IF (Nknots < 1) THEN
				PRINT *, "### MSPL_ERROR ### invoked by new_MSPL"
				PRINT *, "Nknots must be greater than zero"
				STOP
			END IF
			IF (m < 0) THEN
				PRINT *, "### MSPL_ERROR ### invoked by new_MSPL"
				PRINT *, "m must be equal or greater than zero"
				STOP
			END IF
         IF (La>=Lb) THEN
            PRINT *, "### MSPL_ERROR ### invoked by new_MSPL"
            PRINT *, "LB must be greater than LA"
            STOP
         END IF
		END IF
		
		spl%Nknots=Nknots
		
		spl%m=m
		spl%po=2*m+1

		spl%La=La
		spl%Lb=Lb
		spl%delta=(Lb-La)/REAL(Nknots,8)
		spl%Idelta=1.d0/spl%delta

		ALLOCATE(spl%x(-1:Nknots+1))
		DO i = -1, spl%Nknots+1, 1
			spl%x(i)=La+spl%delta*REAL(i,8)
		END DO

		ALLOCATE(spl%t(0:m,0:Nknots))
		spl%t=0.d0
		
		ALLOCATE(spl%S(0:m,0:spl%po))
		CALL generate_matrix_S(spl)
		ALLOCATE(spl%ST(0:spl%po,0:m))
		DO j = 0, m, 1
			DO i = 0, spl%po, 1
				spl%ST(i,j)=spl%S(j,i)
			END DO
		END DO

		ALLOCATE(spl%DF(0:spl%po,0:spl%po))
		CALL generate_derivative_factors(spl)
		
      spl%cutoff=cutoff
      IF (cutoff) THEN
         spl%fade=1
      ELSE
         spl%fade=0
      END IF

		spl%flag_init=.TRUE.
	
	END SUBROUTINE MSPL_new


	SUBROUTINE generate_matrix_S(spl)
		IMPLICIT NONE
		TYPE(MSPLINE) :: spl
		INTEGER :: i1, i2, j, info
		REAL(KIND=8) :: M(0:spl%m,0:spl%m)
		REAL(KIND=8) :: pvt(0:spl%m), work(0:3*spl%m+2)
		
		DO i2 = 0, spl%m, 1
			DO i1 = 0, spl%m, 1
				M(i1,i2)=fact(spl%m+1+i1)/fact(spl%m+1+i1-i2)
			END DO
		END DO
		CALL DGETRF( spl%m+1, spl%m+1, M, spl%m+1, pvt, info )
		IF (info/=0) THEN
			PRINT *, "### MSPL_ERROR ### invoked by generate_matrix_S"
			PRINT *, "Inversion of the matrix M did not work (part1)"
			STOP
		END IF
		CALL DGETRI( spl%m+1, M, spl%m+1, pvt, work, 3*spl%m+3, info )
		IF (info/=0) THEN
			PRINT *, "### MSPL_ERROR ### invoked by generate_matrix_S"
			PRINT *, "Inversion of the matrix M did not work (part2)"
			STOP
		END IF
		
		spl%S=0.d0
		DO i2 = 0, spl%m, 1
			spl%S(i2,i2)=1.d0/(REAL(fact(i2),8))
		END DO
		DO i2 = 0, spl%m, 1
			DO i1 = 0, spl%m, 1
				DO j = 0, i1, 1
					spl%S(i1,i2+spl%m+1)=spl%S(i1,i2+spl%m+1)-M(j,i2)/REAL(fact(i1-j),8)
				END DO
			END DO
		END DO

	END SUBROUTINE generate_matrix_S


	RECURSIVE FUNCTION fact(n) RESULT(f)
		IMPLICIT NONE
		INTEGER(KIND=8) :: f
		INTEGER :: n
		
		IF (n==0) THEN
			f=1
		ELSE
			f=fact(n-1)*n
		END IF
	
	END FUNCTION fact


	SUBROUTINE generate_derivative_factors(spl)
		IMPLICIT NONE
		TYPE(MSPLINE) :: spl
		INTEGER :: der, n, n2, i

		spl%DF(0:spl%po,0)=1.d0
		DO der = 1, spl%po, 1
			DO n = 0, spl%po, 1
				IF (der>n) THEN
					spl%DF(n,der)=0.d0
				ELSE
					spl%DF(n,der)=(fact(n)/fact(n-der))*((-1.d0/spl%delta)**der)
				END IF
			END DO
		END DO

	END SUBROUTINE generate_derivative_factors


	SUBROUTINE MSPL_fit_function(spl,f)
	!generate the spline's parameters which fit the provided function and its first m derivatives
	!the function is supposed to provide REAL(KIND=8) output and to be one-dimensional
		IMPLICIT NONE
		TYPE(MSPLINE) :: spl    !spline to set
		!f(0,x)=function of x   f(1,x)=first derivative computed in x   ...etc
		INTERFACE   !explicit interface for the function f
			FUNCTION f(i,x)
				REAL(KIND=8) :: f
				INTEGER, INTENT(IN) :: i
				REAL(KIND=8), INTENT(IN) :: x
			END FUNCTION f
		END INTERFACE
		INTEGER :: i, j
		
		DO j = 0, spl%Nknots, 1
			DO i = 0, spl%m, 1
				spl%t(i,j)=f(i,spl%x(j))
			END DO
		END DO

	END SUBROUTINE MSPL_fit_function


	FUNCTION compute_MSPL(spl,deriv,r)
	!compute the value in r of the cubic spline specified, and its first deriv derivatives
		IMPLICIT NONE
		REAL(KIND=8) :: compute_MSPL
		INTEGER,INTENT(IN) :: deriv    !specifies the order of the derivative which should be computed
		TYPE(MSPLINE), INTENT(IN) :: spl   !spline to be computed
		REAL(KIND=8), INTENT(IN) :: r   !positione where the spline should be evaluated
		INTEGER :: n, alpha
		INTEGER :: j1, j2    !the two knots between which r is located
		REAL(KIND=8) :: ratio
		REAL(KIND=8) :: r_min_ri, ri_min_r   !r_min_ri=(r-ri)/delta      ri_min_r=(ri-r)/delta
		REAL(KIND=8) :: cia
		REAL(KIND=8) :: R_deriv, dummy1, Rn, R_fact(0:spl%po)          !factors introduced for performance
		REAL(KIND=8) :: t1, t2

		IF (MSPL_DEBUG_MODE) THEN
			IF (.NOT. spl%flag_init) THEN
				PRINT *, "### MSPL_ERROR ###  invoked by compute_MSPL"
				PRINT *, "The spline has not been initialized"
				STOP
			END IF
			IF ( r < spl%La-spl%delta ) THEN
				PRINT *, "### MSPL_ERROR ###  invoked by compute_MSPL"
				PRINT *, "R is too small, outside the range"
				STOP
			END IF
			IF ( r > spl%Lb+spl%delta ) THEN
				PRINT *, "### MSPL_ERROR ###  invoked by compute_MSPL"
				PRINT *, "R is too large, outside the range"
				STOP
			END IF
         IF (deriv<0) THEN
            PRINT *, "### MSPL_ERROR ###  invoked by compute_MSPL"
            PRINT *, "DERIV must be greater than zero"
            STOP
         END IF
         IF (deriv>spl%po) THEN
            PRINT *, "### MSPL_ERROR ###  invoked by compute_MSPL"
            PRINT *, "DERIV is larger than the polynomial order of the spline"
         END IF
		END IF
		
		ratio=r/spl%delta
		j1=FLOOR(ratio)
      !CALL find_index(spl%Nknots,spl%x(0:spl%Nknots),r,j1)
		j2=j1+1
		r_min_ri=(r-spl%x(j1))*spl%Idelta
		ri_min_r=(spl%x(j2)-r)*spl%Idelta

		compute_MSPL=0.d0
		IF ( j1 <= spl%Nknots-spl%fade ) THEN !the -1 make the function fade to 0 when it gets to Lb
			DO alpha = 0, spl%m, 1
				cia=0.d0
				IF (alpha-deriv >= 0) cia=cia+spl%ST(alpha,alpha)*((r_min_ri)**(alpha-deriv))*DABS(spl%DF(alpha,deriv))
				DO n = spl%m+1, spl%po, 1
					IF (n-deriv >= 0) cia=cia+spl%ST(n,alpha)*((r_min_ri)**(n-deriv))*DABS(spl%DF(n,deriv))
				END DO
				compute_MSPL=compute_MSPL+spl%t(alpha,j1)*cia*(spl%delta**(alpha))
			END DO
		END IF
		IF ( j2 <= spl%Nknots-spl%fade ) THEN !the -1 make the function fade to 0 when it gets to Lb
			DO alpha = 0, spl%m, 1
				cia=0.d0
				IF (alpha-deriv >= 0) cia=cia+spl%ST(alpha,alpha)*((ri_min_r)**(alpha-deriv))*spl%DF(alpha,deriv)
				DO n = spl%m+1, spl%po, 1
					IF (n-deriv >= 0) cia=cia+spl%ST(n,alpha)*((ri_min_r)**(n-deriv))*spl%DF(n,deriv)
				END DO
				compute_MSPL=compute_MSPL+spl%t(alpha,j2)*cia*((-spl%delta)**(alpha))
			END DO
		END IF

	END FUNCTION compute_MSPL

   !test done for speed-up, avoiding a recasting from REAL to INTEGER
   !however, there is no gain in performance
   !SUBROUTINE find_index(n,x,r,i)
   !   IMPLICIT NONE
   !   LOGICAL :: flag
   !   INTEGER :: i, i1, i2, n
   !   REAL(KIND=8) :: r, x(0:n)

   !   i1=0
   !   i2=n
   !   i=n/2
   !   flag=.TRUE.
   !   IF (i2-i1<=1) THEN
   !      flag=.FALSE.
   !   END IF
   !   DO WHILE (flag)
   !      IF (r<x(i)) THEN
   !         i2=i
   !      ELSE
   !         i1=i
   !      END IF
   !      i=(i1+i2)/2
   !      IF (i2-i1<=1) THEN
   !         flag=.FALSE.
   !      END IF
   !   END DO

   !END SUBROUTINE find_index
	
	SUBROUTINE MSPL_compute(spl,deriv,r,val,reset)
	!compute the value in r of the cubic spline specified, and its first deriv derivatives
		IMPLICIT NONE
		REAL(KIND=8), INTENT(INOUT) :: val
      LOGICAL, OPTIONAL, INTENT(IN) :: reset
		INTEGER,INTENT(IN) :: deriv    !specifies the order of the derivative which should be computed
		TYPE(MSPLINE), INTENT(IN) :: spl   !spline to be computed
		REAL(KIND=8), INTENT(IN) :: r   !positione where the spline should be evaluated
		INTEGER :: n, alpha
		INTEGER :: j1, j2    !the two knots between which r is located
		REAL(KIND=8) :: ratio
		REAL(KIND=8) :: r_min_ri, ri_min_r   !r_min_ri=(r-ri)/delta      ri_min_r=(ri-r)/delta
		REAL(KIND=8) :: cia
		REAL(KIND=8) :: R_deriv, dummy1, Rn, R_fact(0:spl%po)          !factors introduced for performance
		REAL(KIND=8) :: t1, t2

		IF (MSPL_DEBUG_MODE) THEN
			IF (.NOT. spl%flag_init) THEN
				PRINT *, "### MSPL_ERROR ###  invoked by MSPL_compute"
				PRINT *, "The spline has not been initialized"
				STOP
			END IF
         IF (deriv<0) THEN
            PRINT *, "### MSPL_ERROR ###  invoked by MSPL_compute"
            PRINT *, "DERIV must be greater than zero"
            STOP
         END IF
         IF (deriv>spl%po) THEN
            PRINT *, "### MSPL_ERROR ### invoked by MSPL_compute"
            PRINT *, "DERIV is larger than the polynomial order of the spline"
         END IF
		END IF
		
		ratio=r/spl%delta
		j1=FLOOR(ratio)
		j2=j1+1
		r_min_ri=(r-spl%x(j1))*spl%Idelta
		ri_min_r=(spl%x(j2)-r)*spl%Idelta

      IF (PRESENT(reset)) THEN
         IF (reset) val=0.d0
      ELSE
         val=0.d0
      END IF
		
      IF ( j1 <= spl%Nknots-spl%fade ) THEN !the -1 make the function fade to 0 when it gets to Lb
			DO alpha = 0, spl%m, 1
				cia=0.d0
				IF (alpha-deriv >= 0) cia=cia+spl%ST(alpha,alpha)*((r_min_ri)**(alpha-deriv))*DABS(spl%DF(alpha,deriv))
				!first development, slightly slower
            !DO n = spl%m+1, spl%po, 1
				!	IF (n-deriv >= 0) cia=cia+spl%ST(n,alpha)*((r_min_ri)**(n-deriv))*DABS(spl%DF(n,deriv))
				DO n = MAX(spl%m+1,0), spl%po, 1
					cia=cia+spl%ST(n,alpha)*((r_min_ri)**(n-deriv))*DABS(spl%DF(n,deriv))
				END DO
				val=val+spl%t(alpha,j1)*cia*(spl%delta**(alpha))
			END DO
		END IF
		IF ( j2 <= spl%Nknots-spl%fade ) THEN !the -1 make the function fade to 0 when it gets to Lb
			DO alpha = 0, spl%m, 1
				cia=0.d0
				IF (alpha-deriv >= 0) cia=cia+spl%ST(alpha,alpha)*((ri_min_r)**(alpha-deriv))*spl%DF(alpha,deriv)
				!DO n = spl%m+1, spl%po, 1
				!	IF (n-deriv >= 0) cia=cia+spl%ST(n,alpha)*((ri_min_r)**(n-deriv))*spl%DF(n,deriv)
            DO n = MAX(spl%m+1,0), spl%po, 1
					cia=cia+spl%ST(n,alpha)*((ri_min_r)**(n-deriv))*spl%DF(n,deriv)
				END DO
				val=val+spl%t(alpha,j2)*cia*((-spl%delta)**(alpha))
			END DO
		END IF

	END SUBROUTINE MSPL_compute


   SUBROUTINE MSPL_print_on_file(spl,deriv,filename,npoints)
      IMPLICIT NONE
      TYPE(MSPLINE), INTENT(IN) :: spl
      CHARACTER (LEN=*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: deriv, npoints
      INTEGER :: i
      REAL(KIND=8) :: delta, x, val

		IF (MSPL_DEBUG_MODE) THEN
			IF (deriv < 0) THEN
				PRINT *, "### MSPL_ERROR ### invoked by MSPL_print_on_file"
				PRINT *, "deriv must be greater than zero"
				STOP
			END IF
         IF (npoints<2) THEN
            PRINT *, "### MSPL_ERROR ### invoked by MSPL_print_on_file"
            PRINT *, "npoints must be greater or equal to two"
            STOP
         END IF
		END IF

      OPEN(UNIT=159,FILE=filename,STATUS='UNKNOWN',POSITION='ASIS',IOSTAT=i)
      IF (i/=0) THEN
				PRINT *, "### MSPL_ERROR ### invoked by MSPL_print_on_file"
				PRINT *, "error in opening the file"
				STOP
      END IF
      
      delta=(spl%Lb-spl%La)/REAL(npoints-1,8)
      DO i = 0, npoints-1, 1
         x=i*delta
         CALL MSPL_compute(spl,deriv,x,val)
         WRITE(UNIT=159, FMT=*), x, val
      END DO

      CLOSE(UNIT=159,IOSTAT=i)
      IF (i/=0) THEN
				PRINT *, "### MSPL_ERROR ### invoked by MSPL_print_on_file"
				PRINT *, "error in closing the file"
				STOP
      END IF
      
   END SUBROUTINE MSPL_print_on_file


	SUBROUTINE MSPL_t_deriv(spl,r,t_deriv,reset)
	!compute the m*(Nknots+1) variational derivatives
		IMPLICIT NONE
		TYPE(MSPLINE), INTENT(IN) :: spl
		REAL(KIND=8), INTENT(IN) :: r
      LOGICAL, OPTIONAL, INTENT(IN) :: reset           !accumulate the new values in t_deriv?
		INTEGER :: j1, j2, alpha, n
		REAL(KIND=8) :: ratio, r_min_ri, ri_min_r, cia
		REAL(KIND=8), INTENT(INOUT) :: t_deriv(0:spl%m,0:spl%Nknots)

      IF (MSPL_DEBUG_MODE) THEN
		   IF (.NOT. spl%flag_init) THEN
		   	PRINT *, "### MSPL_ERROR ### invoked by MSPL_t_deriv"
		   	PRINT *, "The spline has not been initialized"
		   	STOP
		   END IF
			IF ( r < spl%La-spl%delta ) THEN
				PRINT *, "### MSPL_ERROR ### invoked by MSPL_t_deriv"
				PRINT *, "R is too small, outside the range"
				STOP
			END IF
			IF ( r > spl%Lb+spl%delta ) THEN
				PRINT *, "### MSPL_ERROR ### invoked by MSPL_t_deriv"
				PRINT *, "R is too large, outside the range"
				STOP
			END IF
      END IF
		
		ratio=r/spl%delta
		j1=FLOOR(ratio)
		j2=j1+1
		r_min_ri=(r-spl%x(j1))*spl%Idelta
		ri_min_r=(spl%x(j2)-r)*spl%Idelta

      IF (.NOT. PRESENT(reset)) THEN
         t_deriv=0.d0
      ELSE
         IF (reset) t_deriv=0.d0
      END IF

		IF ( j1 <= spl%Nknots-spl%fade ) THEN
			DO alpha = 0, spl%m, 1
				cia=spl%ST(alpha,alpha)*(r_min_ri**alpha)
				DO n = spl%m+1, spl%po, 1
					cia=cia+spl%ST(n,alpha)*(r_min_ri**n)
				END DO
				t_deriv(alpha,j1)=t_deriv(alpha,j1)+cia*(spl%delta**alpha)
			END DO
		END IF
		IF ( j2 <= spl%Nknots-spl%fade ) THEN
			DO alpha = 0, spl%m, 1
				cia=spl%ST(alpha,alpha)*(ri_min_r**alpha)
				DO n = spl%m+1, spl%po, 1
					cia=cia+spl%ST(n,alpha)*(ri_min_r**n)
				END DO
				t_deriv(alpha,j2)=t_deriv(alpha,j2)+cia*((-spl%delta)**alpha)
			END DO
		END IF
	
      END SUBROUTINE MSPL_t_deriv

	
	SUBROUTINE MSPL_boxNpart_t_deriv(spl,R,ndim,Npart,Lbox,t_deriv,reset)
		IMPLICIT NONE
		TYPE(MSPLINE), INTENT(IN) :: spl
		INTEGER, INTENT(IN) :: Npart, ndim
      LOGICAL, OPTIONAL, INTENT(IN) :: reset
		REAL(KIND=8), INTENT(IN) :: R(1:ndim,1:Npart), Lbox(1:ndim)
		INTEGER :: j1, j2, alpha, n, ii, ij
		REAL(KIND=8) :: ratio, r_min_ri, ri_min_r, cia, dist(0:3)
		REAL(KIND=8), INTENT(OUT) :: t_deriv(0:spl%m,0:spl%Nknots)

      IF (MSPL_DEBUG_MODE) THEN
		   IF (.NOT. spl%flag_init) THEN
		   	PRINT *, "### MSPL_ERROR ### invoked by MSPL_t_Npart_deriv"
		   	PRINT *, "The spline has not been initialized"
		   	STOP
		   END IF
         IF (ndim<1) THEN
            PRINT *, "### MSPL_ERROR ### invoked by MSPL_t_Npart_deriv"
            PRINT *, "NDIM must be greater or equal to zero"
            STOP
         END IF
         IF (Npart<1) THEN
            PRINT *, "### MSPL_ERROR ### invoked by MSPL_t_Npart_deriv"
            PRINT *, "NPART must be greater or equal to zero"
            STOP
         END IF
         IF (MINVAL(Lbox(1:ndim))<=0.d0) THEN
            PRINT *, "### MSPL_ERROR ### invoked by MSPL_t_Npart_deriv"
            PRINT *, "All the components of LBOX must be greater than zero"
            STOP
         END IF
      END IF

      IF (.NOT. PRESENT(reset)) THEN
         t_deriv=0.d0
      ELSE
         IF (reset) t_deriv=0.d0
      END IF

		DO ii = 1, Npart-1, 1
			DO ij = ii+1, Npart, 1
				dist(1:ndim)=R(1:ndim,ij)-R(1:ndim,ii)
            dist(1:ndim)=dist(1:ndim)-Lbox(1:ndim)*DNINT(dist(1:ndim)/Lbox(1:ndim))
				dist(0)=DSQRT(DOT_PRODUCT(dist(1:ndim),dist(1:ndim)))
            IF (dist(0) < spl%Lb) THEN
				   ratio=dist(0)/spl%delta
				   j1=FLOOR(ratio)
				   j2=j1+1
				   r_min_ri=(dist(0)-spl%x(j1))*spl%Idelta
				   ri_min_r=(spl%x(j2)-dist(0))*spl%Idelta

				   IF ( j1 <= spl%Nknots-spl%fade ) THEN
				   	DO alpha = 0, spl%m, 1
				   		cia=spl%ST(alpha,alpha)*(r_min_ri**alpha)
				   		DO n = spl%m+1, spl%po, 1
				   			cia=cia+spl%ST(n,alpha)*(r_min_ri**n)
				   		END DO
				   		t_deriv(alpha,j1)=t_deriv(alpha,j1)+cia*(spl%delta**alpha)
				   	END DO
				   END IF
				   IF ( j2 <= spl%Nknots-spl%fade ) THEN
				   	DO alpha = 0, spl%m, 1
				   		cia=spl%ST(alpha,alpha)*(ri_min_r**alpha)
				   		DO n = spl%m+1, spl%po, 1
				   			cia=cia+spl%ST(n,alpha)*(ri_min_r**n)
				   		END DO
				   		t_deriv(alpha,j2)=t_deriv(alpha,j2)+cia*((-spl%delta)**alpha)
				   	END DO
				   END IF
            END IF
			END DO
		END DO
		
	END SUBROUTINE MSPL_boxNpart_t_deriv


   SUBROUTINE MSPL_grad_t_deriv(spl,r_vec,ndim,grad_t_deriv,reset)
      IMPLICIT NONE
      TYPE(MSPLINE), INTENT(IN) :: spl
      INTEGER, INTENT(IN) :: ndim
      REAL(KIND=8), INTENT(IN) :: r_vec(0:ndim)
      LOGICAL, OPTIONAL, INTENT(IN) :: reset
      REAL(KIND=8), INTENT(INOUT) :: grad_t_deriv(0:spl%m,0:spl%Nknots,1:ndim)
      INTEGER :: j1, j2, alpha, n
      REAL(KIND=8) :: r, ratio, r_min_ri, ri_min_r, cia, foo
      REAL(KIND=8) :: grad_r(0:spl%m,1:2)

      IF (MSPL_DEBUG_MODE) THEN
		   IF (.NOT. spl%flag_init) THEN
		   	PRINT *, "### MSPL_ERROR ### invoked by MSPL_t_Npart_deriv"
		   	PRINT *, "The spline has not been initialized"
		   	STOP
		   END IF
         IF (ndim<1) THEN
            PRINT *, "### MSPL_ERROR ### invoked by MSPL_t_Npart_deriv"
            PRINT *, "NDIM must be greater or equal to zero"
            STOP
         END IF
      END IF

		ratio=r_vec(0)/spl%delta
		j1=FLOOR(ratio)
		j2=j1+1
		r_min_ri=(r_vec(0)-spl%x(j1))*spl%Idelta
		ri_min_r=(spl%x(j2)-r_vec(0))*spl%Idelta

      !compute the radial gradient
      grad_r=0.d0
		IF ( j1 <= spl%Nknots-spl%fade ) THEN
			DO alpha = 0, spl%m, 1
				cia=0.d0
            IF (alpha-1 >= 0) cia=cia+spl%ST(alpha,alpha)*((r_min_ri)**(alpha-1))*DABS(spl%DF(alpha,1))
				DO n = MAX(spl%m+1,0), spl%po, 1
               cia=cia+spl%ST(n,alpha)*((r_min_ri)**(n-1))*DABS(spl%DF(n,1))
				END DO
            grad_r(alpha,1)=grad_r(alpha,1)+cia*(spl%delta**alpha)
			END DO
		END IF
		IF ( j2 <= spl%Nknots-spl%fade ) THEN
			DO alpha = 0, spl%m, 1
            cia=0.d0
            IF (alpha-1 >=0) cia=cia+spl%ST(alpha,alpha)*((ri_min_r)**(alpha-1))*spl%DF(alpha,1)
				DO n = MAX(spl%m+1,0), spl%po, 1
               cia=cia+spl%ST(n,alpha)*((ri_min_r)**(n-1))*spl%DF(n,1)
				END DO
            grad_r(alpha,2)=grad_r(alpha,2)+cia*((-spl%delta)**alpha)
			END DO
		END IF

      !reset the gradient, if not asked to do otherwise
      IF (PRESENT(reset)) THEN
         IF (reset) grad_t_deriv=0.d0
      ELSE
         grad_t_deriv=0.d0
      END IF

      !update the gradient, adding the vectorial component
      DO n = 1, ndim, 1
         grad_t_deriv(:,j1:j2,n)=grad_t_deriv(:,j1:j2,n)+grad_r(:,:)*r_vec(n)/r_vec(0)
      END DO

   END SUBROUTINE MSPL_grad_t_deriv


   SUBROUTINE MSPL_boxNpart_grad_t_deriv(spl,R,ndim,npart,Lbox,grad_t_deriv,reset)
      IMPLICIT NONE
      TYPE(MSPLINE), INTENT(IN) :: spl
      INTEGER, INTENT(IN) :: ndim, npart
      REAL(KIND=8), INTENT(IN) :: R(1:ndim,1:npart), Lbox(1:ndim)
      REAL(KIND=8), INTENT(INOUT) :: grad_t_deriv(0:spl%m,0:spl%Nknots,1:ndim,1:npart)
      LOGICAL, OPTIONAL, INTENT(IN) :: reset
      INTEGER :: i1, i2, j1,j2,j3
      REAL(KIND=8) :: dist(0:ndim)
      REAL(KIND=8) :: grad(0:spl%m,0:spl%Nknots,1:ndim)

      IF (MSPL_DEBUG_MODE) THEN
         IF (Npart<1) THEN
            PRINT *, "### MSPL_ERROR ### invoked by MSPL_t_Npart_deriv"
            PRINT *, "NPART must be greater or equal to zero"
            STOP
         END IF
         IF (MINVAL(Lbox(1:ndim))<=0.d0) THEN
            PRINT *, "### MSPL_ERROR ### invoked by MSPL_t_Npart_deriv"
            PRINT *, "All the components of LBOX must be greater than zero"
            STOP
         END IF
      END IF

      IF (PRESENT(reset)) THEN
         IF (reset) grad_t_deriv=0.d0
      ELSE
        grad_t_deriv=0.d0 
      END IF

      DO i2 = 1, npart-1, 1
      DO i1 = i2+1, npart, 1
			dist(1:ndim)=R(1:ndim,i2)-R(1:ndim,i1)
         dist(1:ndim)=dist(1:ndim)-Lbox(1:ndim)*DNINT(dist(1:ndim)/Lbox(1:ndim))
			dist(0)=DSQRT(DOT_PRODUCT(dist(1:ndim),dist(1:ndim)))
         CALL MSPL_grad_t_deriv(spl,dist(0:ndim),ndim,grad(:,:,:))
         grad_t_deriv(:,:,:,i2)=grad_t_deriv(:,:,:,i2)+grad(:,:,:)
         grad_t_deriv(:,:,:,i1)=grad_t_deriv(:,:,:,i1)-grad(:,:,:)
      END DO
      END DO
      
   END SUBROUTINE MSPL_boxNpart_grad_t_deriv


   SUBROUTINE MSPL_lapl_t_deriv(spl,r_vec,ndim,lapl_t_deriv,reset)
      IMPLICIT NONE
      TYPE(MSPLINE), INTENT(IN) :: spl
      INTEGER, INTENT(IN) :: ndim
      REAL(KIND=8), INTENT(IN) :: r_vec(0:ndim)
      LOGICAL, OPTIONAL, INTENT(IN) :: reset
      REAL(KIND=8), INTENT(INOUT) :: lapl_t_deriv(0:spl%m,0:spl%Nknots)
      INTEGER :: j1, j2, alpha, n
      REAL(KIND=8) :: r, ratio, r_min_ri, ri_min_r, cia, foo
      REAL(KIND=8) :: grad_r(0:spl%m,1:2), lapl_r(0:spl%m,1:2)

      IF (MSPL_DEBUG_MODE) THEN
		   IF (.NOT. spl%flag_init) THEN
		   	PRINT *, "### MSPL_ERROR ### invoked by MSPL_t_Npart_deriv"
		   	PRINT *, "The spline has not been initialized"
		   	STOP
		   END IF
         IF (ndim<1) THEN
            PRINT *, "### MSPL_ERROR ### invoked by MSPL_t_Npart_deriv"
            PRINT *, "NDIM must be greater or equal to zero"
            STOP
         END IF
      END IF

		ratio=r_vec(0)/spl%delta
		j1=FLOOR(ratio)
		j2=j1+1
		r_min_ri=(r_vec(0)-spl%x(j1))*spl%Idelta
		ri_min_r=(spl%x(j2)-r_vec(0))*spl%Idelta

      !compute the radial gradient
      grad_r=0.d0
		IF ( j1 <= spl%Nknots-spl%fade ) THEN
			DO alpha = 0, spl%m, 1
				cia=0.d0
            IF (alpha-1 >= 0) cia=cia+spl%ST(alpha,alpha)*((r_min_ri)**(alpha-1))*DABS(spl%DF(alpha,1))
				DO n = MAX(spl%m+1,0), spl%po, 1
               cia=cia+spl%ST(n,alpha)*((r_min_ri)**(n-1))*DABS(spl%DF(n,1))
				END DO
            grad_r(alpha,1)=grad_r(alpha,1)+cia*(spl%delta**alpha)
			END DO
		END IF
		IF ( j2 <= spl%Nknots-spl%fade ) THEN
			DO alpha = 0, spl%m, 1
            cia=0.d0
            IF (alpha-1 >=0) cia=cia+spl%ST(alpha,alpha)*((ri_min_r)**(alpha-1))*spl%DF(alpha,1)
				DO n = MAX(spl%m+1,0), spl%po, 1
               cia=cia+spl%ST(n,alpha)*((ri_min_r)**(n-1))*spl%DF(n,1)
				END DO
            grad_r(alpha,2)=grad_r(alpha,2)+cia*((-spl%delta)**alpha)
			END DO
		END IF
      !compute the radial laplacian
      lapl_r=0.d0
		IF ( j1 <= spl%Nknots-spl%fade ) THEN
			DO alpha = 0, spl%m, 1
				cia=0.d0
            IF (alpha-2 >= 0) cia=cia+spl%ST(alpha,alpha)*((r_min_ri)**(alpha-2))*DABS(spl%DF(alpha,2))
				DO n = MAX(spl%m+1,0), spl%po, 1
               cia=cia+spl%ST(n,alpha)*((r_min_ri)**(n-2))*DABS(spl%DF(n,2))
				END DO
            lapl_r(alpha,1)=lapl_r(alpha,1)+cia*(spl%delta**alpha)
			END DO
		END IF
		IF ( j2 <= spl%Nknots-spl%fade ) THEN
			DO alpha = 0, spl%m, 1
            cia=0.d0
            IF (alpha-2 >=0) cia=cia+spl%ST(alpha,alpha)*((ri_min_r)**(alpha-2))*spl%DF(alpha,2)
				DO n = MAX(spl%m+1,0), spl%po, 1
               cia=cia+spl%ST(n,alpha)*((ri_min_r)**(n-2))*spl%DF(n,2)
				END DO
            lapl_r(alpha,2)=lapl_r(alpha,2)+cia*((-spl%delta)**alpha)
			END DO
		END IF

      !reset the laplacian, if not asked to do otherwise
      IF (PRESENT(reset)) THEN
         IF (reset) lapl_t_deriv=0.d0
      ELSE
         lapl_t_deriv=0.d0
      END IF

      !update the laplacian, adding the vectorial component
      ratio=1.d0/r_vec(0)
      DO n = 1, ndim, 1
         lapl_t_deriv(:,j1:j2)=lapl_t_deriv(:,j1:j2)+ratio*( grad_r(:,1:2) + &
                                 ratio*r_vec(n)*r_vec(n)*( lapl_r(:,1:2) - &
                                 ratio*grad_r(:,1:2) ) )
      END DO

   END SUBROUTINE MSPL_lapl_t_deriv


   SUBROUTINE MSPL_grad_and_lapl_t_deriv(spl,r_vec,ndim,grad_t_deriv,lapl_t_deriv,reset)
      IMPLICIT NONE
      TYPE(MSPLINE), INTENT(IN) :: spl
      INTEGER, INTENT(IN) :: ndim
      REAL(KIND=8), INTENT(IN) :: r_vec(0:ndim)
      LOGICAL, OPTIONAL, INTENT(IN) :: reset
      REAL(KIND=8), INTENT(INOUT) :: grad_t_deriv(0:spl%m,0:spl%Nknots,1:ndim)
      REAL(KIND=8), INTENT(INOUT) :: lapl_t_deriv(0:spl%m,0:spl%Nknots) 
      INTEGER :: j1, j2, alpha, n
      REAL(KIND=8) :: r, ratio, r_min_ri, ri_min_r, cia, foo
      REAL(KIND=8) :: grad_r(0:spl%m,1:2), lapl_r(0:spl%m,1:2)

      IF (MSPL_DEBUG_MODE) THEN
		   IF (.NOT. spl%flag_init) THEN
		   	PRINT *, "### MSPL_ERROR ### invoked by MSPL_t_Npart_deriv"
		   	PRINT *, "The spline has not been initialized"
		   	STOP
		   END IF
         IF (ndim<1) THEN
            PRINT *, "### MSPL_ERROR ### invoked by MSPL_t_Npart_deriv"
            PRINT *, "NDIM must be greater or equal to zero"
            STOP
         END IF
      END IF

		ratio=r_vec(0)/spl%delta
		j1=FLOOR(ratio)
		j2=j1+1
		r_min_ri=(r_vec(0)-spl%x(j1))*spl%Idelta
		ri_min_r=(spl%x(j2)-r_vec(0))*spl%Idelta

      !compute the radial gradient
      grad_r=0.d0
		IF ( j1 <= spl%Nknots-spl%fade ) THEN
			DO alpha = 0, spl%m, 1
				cia=0.d0
            IF (alpha-1 >= 0) cia=cia+spl%ST(alpha,alpha)*((r_min_ri)**(alpha-1))*DABS(spl%DF(alpha,1))
				DO n = MAX(spl%m+1,0), spl%po, 1
               cia=cia+spl%ST(n,alpha)*((r_min_ri)**(n-1))*DABS(spl%DF(n,1))
				END DO
            grad_r(alpha,1)=grad_r(alpha,1)+cia*(spl%delta**alpha)
			END DO
		END IF
		IF ( j2 <= spl%Nknots-spl%fade ) THEN
			DO alpha = 0, spl%m, 1
            cia=0.d0
            IF (alpha-1 >=0) cia=cia+spl%ST(alpha,alpha)*((ri_min_r)**(alpha-1))*spl%DF(alpha,1)
				DO n = MAX(spl%m+1,0), spl%po, 1
               cia=cia+spl%ST(n,alpha)*((ri_min_r)**(n-1))*spl%DF(n,1)
				END DO
            grad_r(alpha,2)=grad_r(alpha,2)+cia*((-spl%delta)**alpha)
			END DO
		END IF
      !compute the radial laplacian
      lapl_r=0.d0
		IF ( j1 <= spl%Nknots-spl%fade ) THEN
			DO alpha = 0, spl%m, 1
				cia=0.d0
            IF (alpha-2 >= 0) cia=cia+spl%ST(alpha,alpha)*((r_min_ri)**(alpha-2))*DABS(spl%DF(alpha,2))
				DO n = MAX(spl%m+1,0), spl%po, 1
               cia=cia+spl%ST(n,alpha)*((r_min_ri)**(n-2))*DABS(spl%DF(n,2))
				END DO
            lapl_r(alpha,1)=lapl_r(alpha,1)+cia*(spl%delta**alpha)
			END DO
		END IF
		IF ( j2 <= spl%Nknots-spl%fade ) THEN
			DO alpha = 0, spl%m, 1
            cia=0.d0
            IF (alpha-2 >=0) cia=cia+spl%ST(alpha,alpha)*((ri_min_r)**(alpha-2))*spl%DF(alpha,2)
				DO n = MAX(spl%m+1,0), spl%po, 1
               cia=cia+spl%ST(n,alpha)*((ri_min_r)**(n-2))*spl%DF(n,2)
				END DO
            lapl_r(alpha,2)=lapl_r(alpha,2)+cia*((-spl%delta)**alpha)
			END DO
		END IF

      !reset the gradient, if not asked to do otherwise
      IF (PRESENT(reset)) THEN
         IF (reset) THEN
            grad_t_deriv=0.d0
         END IF
      ELSE
         grad_t_deriv=0.d0
      END IF

      !store this value to gain in performance
      ratio=1.d0/r_vec(0)

      !update the gradient
      DO n = 1, ndim, 1
         grad_t_deriv(:,j1:j2,n)=grad_t_deriv(:,j1:j2,n)+grad_r(:,:)*r_vec(n)*ratio
      END DO

      !reset the laplacian, if not asked to do otherwise
      IF (PRESENT(reset)) THEN
         IF (reset) THEN
            lapl_t_deriv=0.d0
         END IF
      ELSE
         lapl_t_deriv=0.d0
      END IF

      !update the laplacian
      DO n = 1, ndim, 1
         lapl_t_deriv(:,j1:j2)=lapl_t_deriv(:,j1:j2)+ratio*( grad_r(:,1:2) + &
                                 ratio*r_vec(n)*r_vec(n)*( lapl_r(:,1:2) - &
                                 ratio*grad_r(:,1:2) ) )
      END DO
      
   END SUBROUTINE MSPL_grad_and_lapl_t_deriv


   SUBROUTINE MSPL_boxNpart_grad_and_lapl_t_deriv(spl,R,ndim,npart,Lbox,grad_t_deriv,lapl_t_deriv,reset)
      IMPLICIT NONE
      TYPE(MSPLINE), INTENT(IN) :: spl
      INTEGER, INTENT(IN) :: ndim, npart
      REAL(KIND=8), INTENT(IN) :: R(1:ndim,1:npart), Lbox(1:ndim)
      REAL(KIND=8), INTENT(INOUT) :: grad_t_deriv(0:spl%m,0:spl%Nknots,1:ndim,1:npart)
      REAL(KIND=8), INTENT(INOUT) :: lapl_t_deriv(0:spl%m,0:spl%Nknots,1:npart)
      LOGICAL, OPTIONAL, INTENT(IN) :: reset 
      INTEGER :: i1, i2
      REAL(KIND=8) :: dist(0:ndim)
      REAL(KIND=8) :: grad(0:spl%m,0:spl%Nknots,1:ndim)
      REAL(KIND=8) :: lapl(0:spl%m,0:spl%Nknots)
      
      IF (MSPL_DEBUG_MODE) THEN
         IF (Npart<1) THEN
            PRINT *, "### MSPL_ERROR ### invoked by MSPL_t_Npart_deriv"
            PRINT *, "NPART must be greater or equal to zero"
            STOP
         END IF
         IF (MINVAL(Lbox(1:ndim))<=0.d0) THEN
            PRINT *, "### MSPL_ERROR ### invoked by MSPL_t_Npart_deriv"
            PRINT *, "All the components of LBOX must be greater than zero"
            STOP
         END IF
      END IF

      IF (PRESENT(reset)) THEN
         IF (reset) THEN
            grad_t_deriv=0.d0
            lapl_t_deriv=0.d0
         END IF
      ELSE
        grad_t_deriv=0.d0 
        lapl_t_deriv=0.d0
      END IF

      DO i2 = 1, npart-1, 1
      DO i1 = i2+1, npart, 1
			dist(1:ndim)=R(1:ndim,i2)-R(1:ndim,i1)
         dist(1:ndim)=dist(1:ndim)-Lbox(1:ndim)*DNINT(dist(1:ndim)/Lbox(1:ndim))
			dist(0)=DSQRT(DOT_PRODUCT(dist(1:ndim),dist(1:ndim)))
         CALL MSPL_grad_and_lapl_t_deriv(spl,dist,ndim,grad,lapl)
         grad_t_deriv(:,:,:,i2)=grad_t_deriv(:,:,:,i2)+grad(:,:,:)
         grad_t_deriv(:,:,:,i1)=grad_t_deriv(:,:,:,i1)-grad(:,:,:)
         lapl_t_deriv(:,:,i2)=lapl_t_deriv(:,:,i2)+lapl(:,:)
         lapl_t_deriv(:,:,i1)=lapl_t_deriv(:,:,i1)+lapl(:,:)
      END DO
      END DO
      
   END SUBROUTINE MSPL_boxNpart_grad_and_lapl_t_deriv


	SUBROUTINE MSPL_carbon_copy(original_spl,cc_spl)
		IMPLICIT NONE
		TYPE(MSPLINE), INTENT(IN) :: original_spl
		TYPE(MSPLINE), INTENT(INOUT) :: cc_spl
		INTEGER :: ix, i1
		
      IF (MSPL_DEBUG_MODE) THEN
		   IF (.NOT. original_spl%flag_init) THEN
		   	PRINT *, "### MSPL_ERROR ### invoked by MSPL_carbon_copy"
		   	PRINT *, "The original spline has not been initialized"
		   	STOP
         END IF
		   IF (.NOT. cc_spl%flag_init) THEN
		   	PRINT *, "### MSPL_ERROR ### invoked by MSPL_carbon_copy"
		   	PRINT *, "The carbon-copy spline has not been initialized"
		   	STOP
         END IF
      END IF

		DO ix = 0, cc_spl%Nknots, 1
			DO i1 = 0, cc_spl%m, 1
				cc_spl%t(i1,ix)=compute_MSPL(original_spl,i1,cc_spl%x(ix))
			END DO
		END DO

	END SUBROUTINE MSPL_carbon_copy


	SUBROUTINE MSPL_deallocate(spl)
		IMPLICIT NONE
		TYPE(MSPLINE) :: spl

      IF (MSPL_DEBUG_MODE) THEN
		   IF (.NOT. spl%flag_init) THEN
		   	PRINT *, "### MSPL_ERROR ### invoked by MSPL_deallocate"
		   	PRINT *, "The spline has not been initialized"
		   	STOP
         END IF
      END IF
		
		DEALLOCATE(spl%x)
		DEALLOCATE(spl%t)
		DEALLOCATE(spl%S)
		DEALLOCATE(spl%ST)
		DEALLOCATE(spl%DF)
		
	END SUBROUTINE MSPL_deallocate


   SUBROUTINE MSPL_change_debug_mode(debug_mode)
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: debug_mode

      MSPL_DEBUG_MODE=debug_mode
      
   END SUBROUTINE MSPL_change_debug_mode


   FUNCTION get_debug_mode_MSPL()
      IMPLICIT NONE
      LOGICAL :: get_debug_mode_MSPL

      get_debug_mode_MSPL=MSPL_DEBUG_MODE
   
   END FUNCTION get_debug_mode_MSPL

	
END MODULE markuspline


