INCLUDE 'module_markuspline.f90'
INCLUDE 'module_function.f90'

PROGRAM prova
	USE example
	USE markuspline
	IMPLICIT NONE
	REAL(KIND=8), PARAMETER :: L1=0.d0, L2=10.d0
	INTEGER, PARAMETER :: Nknots1=100, m1=3
	INTEGER, PARAMETER :: Nknots2=5, m2=5
	INTEGER :: i1, i2, NNN
	INTEGER :: j1, j2
	REAL(KIND=8) :: x, y, delta
	REAL(KIND=8) :: t1, t2
	REAL(KIND=8) :: v1, v2, dv
	TYPE(MSPLINE) :: s1, s2
	REAL(KIND=8), ALLOCATABLE :: der_t(:,:), der_t2(:,:)
   INTEGER, PARAMETER :: ndim=3, npart=7
   REAL(KIND=8) :: R(1:ndim,npart), Lbox(1:ndim), dist(0:ndim)
   REAL(KIND=8), ALLOCATABLE :: grad(:,:,:), der_t3(:,:,:)
   REAL(KIND=8), ALLOCATABLE :: der_t4(:,:,:,:), der_t5(:,:,:,:)
   REAL(KIND=8), ALLOCATABLE :: der_t6(:,:), lapl(:,:)
   REAL(KIND=8), ALLOCATABLE :: der_t7(:,:,:)

	s1=new_MSPL(m=m1,La=L1,Nknots=Nknots1,Lb=L2,CUTOFF=.FALSE.)
	CALL MSPL_fit_function(s1,trial_function)
	
	!compute s1 and print on file on a very dense grid
	NNN=1000
   IF (m1>=0) CALL MSPL_print_on_file(s1,0,'vs0.dat',NNN)
   IF (m1>=1) CALL MSPL_print_on_file(s1,1,'vs1.dat',NNN)
   IF (m1>=2) CALL MSPL_print_on_file(s1,2,'vs2.dat',NNN)
   IF (m1>=3) CALL MSPL_print_on_file(s1,3,'vs3.dat',NNN)
   IF (m1>=4) CALL MSPL_print_on_file(s1,4,'vs4.dat',NNN)
   IF (m1>=5) CALL MSPL_print_on_file(s1,5,'vs5.dat',NNN)

	PRINT *, "Check the performance - a lot of calculations!"
	NNN=10
	delta=(L2-L1)/REAL(NNN,8)
	CALL CPU_TIME(t1)
	DO i2 = 1, 10000, 1
		DO i1 = 0, NNN-1, 1
			x=i1*delta
			!y=compute_MSPL(s1,0,x)
         CALL MSPL_compute(s1,2,x,y)
		END DO
	END DO
	CALL CPU_TIME(t2)
	PRINT *, "TIME = ", t2-t1
	PRINT *, 
	
	!create a new MSPLINE s2
	s2=new_MSPL(m2,Nknots2,L1,L2,.TRUE.)

	!make a carbon copy of s1, and save it in s2
	CALL MSPL_carbon_copy(s1,s2)
	
	!CALCULATE s2 on a very dense grid	
	NNN=1000
   IF (m1>=0) CALL MSPL_print_on_file(s2,0,'vss0.dat',NNN)
   IF (m1>=1) CALL MSPL_print_on_file(s2,1,'vss1.dat',NNN)
   IF (m1>=2) CALL MSPL_print_on_file(s2,2,'vss2.dat',NNN)
   IF (m1>=3) CALL MSPL_print_on_file(s2,3,'vss3.dat',NNN)
   IF (m1>=4) CALL MSPL_print_on_file(s2,4,'vss4.dat',NNN)
   IF (m1>=5) CALL MSPL_print_on_file(s2,5,'vss5.dat',NNN)
	
	PRINT *, "> > > Check the variational derivatives (MSPL_t_deriv)"
	!PRINT *, "ranges for t are ", 0, s1%m, " and ", 0, s1%Nknots-1
	dv=0.1d0
	ALLOCATE(der_t(0:s1%m,0:s1%Nknots))
	DO i1 = 1, 3, 1
		CALL RANDOM_NUMBER(x)
		j1=x*(s1%m+1)
		CALL RANDOM_NUMBER(x)
		j2=x*s1%Nknots
		CALL RANDOM_NUMBER(x)
		x=s1%x(j2)+x*(s1%delta)
		!PRINT *, "index of t: ", j1, j2
		!PRINT *, "x = ", x
		v1=compute_MSPL(s1,0,x)
		t1=s1%t(j1,j2)
		s1%t(j1,j2)=s1%t(j1,j2)+dv
		v2=compute_MSPL(s1,0,x)
		s1%t(j1,j2)=t1
		!PRINT *, "v1 and v2 are : ", v1, v2
		PRINT *, "Numerical derivative =  ", (v2-v1)/dv
      CALL MSPL_t_deriv(s1,x,der_t)
		PRINT *, "Analitical derivative = ", der_t(j1,j2)
		PRINT *, 
	END DO

   PRINT *, "> > > Check the all-particles derivative (MSPL_Npart_t_deriv)"
   CALL RANDOM_NUMBER(R)
   R=s1%La+R*(s1%Lb-s1%La)
   Lbox=s1%Lb-s1%La
   ALLOCATE(der_t2(0:s1%m,0:s1%Nknots))
   CALL MSPL_boxNpart_t_deriv(s1,R,ndim,npart,Lbox,der_t2)
   der_t=0.d0
   DO i2 = 1, npart-1, 1
       DO i1 = i2+1, npart, 1
         dist(1:ndim)=R(1:ndim,i2)-R(1:ndim,i1)
         dist(1:ndim)=dist(1:ndim)-Lbox(1:ndim)*DNINT(dist(1:ndim)/Lbox(1:ndim))
         dist(0)=DSQRT(DOT_PRODUCT(dist(1:ndim),dist(1:ndim)))
         !PRINT *, dist(0), Lbox
         CALL MSPL_t_deriv(s1,dist(0),der_t,RESET=.FALSE.)
       END DO
   END DO
   PRINT *, "ALL-PARTICLES: ", SUM(der_t2)
   PRINT *, "ONE-BY-ONE   : ", SUM(der_t)
   PRINT *, "Difference   : ", SUM(der_t2-der_t)
   PRINT *, 

   PRINT *, "> > > Check the gradient of the variational derivative (MSPL_grad_t_deriv)"
   ALLOCATE(der_t3(0:s1%m,0:s1%Nknots,1:ndim))
   ALLOCATE(grad(0:s1%m,0:s1%Nknots,1:ndim))
   dv=0.0000001d0
   DO i2 = 1, 3, 1
		CALL RANDOM_NUMBER(x)
		j1=x*(s1%m+1)
		CALL RANDOM_NUMBER(x)
		j2=x*s1%Nknots
      CALL RANDOM_NUMBER(dist(1:ndim))
      dist(1:ndim)=s1%La+dist(1:ndim)*(s1%Lb-s1%La)*0.5d0
      dist(0)=DSQRT(DOT_PRODUCT(dist(1:ndim),dist(1:ndim)))
      CALL MSPL_t_deriv(s1,dist(0),der_t)
      grad=0.d0
      DO i1 = 1, ndim, 1
         y=dist(i1)
         dist(i1)=dist(i1)+dv
         dist(0)=DSQRT(DOT_PRODUCT(dist(1:ndim),dist(1:ndim)))
         CALL MSPL_t_deriv(s1,dist(0),der_t2)
         grad(:,:,i1)=(der_t2-der_t)/dv
         dist(i1)=y
         dist(0)=DSQRT(DOT_PRODUCT(dist(1:ndim),dist(1:ndim)))
      END DO
      PRINT *, "Numerical gradient = ", SUM(grad(:,:,:))
      grad=0.d0
      CALL MSPL_grad_t_deriv(s1,dist(0:ndim),ndim,grad)
      PRINT *, "Analitical gradient = ", SUM(grad(:,:,:))
      PRINT *, 
   END DO

   PRINT *, "> > > Check the all-particle gradient of the variational derivative (MSPL_boxNpart_grad_t_deriv)"
   CALL RANDOM_NUMBER(R)
   ALLOCATE(der_t4(0:s1%m,0:s1%Nknots,1:ndim,1:npart))
   ALLOCATE(der_t5(0:s1%m,0:s1%Nknots,1:ndim,1:npart))
   der_t4=0.d0
   R=s1%La+R*(s1%Lb-s1%La)
   Lbox=s1%Lb-s1%La
   der_t5=0.d0
   DO i2 = 1, npart, 1
   DO i1 = 1, npart, 1
      IF (i1 /= i2) THEN
         dist(1:ndim)=R(1:ndim,i2)-R(1:ndim,i1)
         dist(1:ndim)=dist(1:ndim)-Lbox(1:ndim)*DNINT(dist(1:ndim)/Lbox(1:ndim))
         dist(0)=DSQRT(DOT_PRODUCT(dist(1:ndim),dist(1:ndim)))
         CALL MSPL_grad_t_deriv(s1,dist(0:ndim),ndim,der_t4(:,:,:,i2))
         der_t5(:,:,:,i2)=der_t5(:,:,:,i2)+der_t4(:,:,:,i2)
         !CALL MSPL_grad_t_deriv(s1,dist(0:ndim),ndim,der_t5(:,:,:,i2),RESET=.FALSE.)
      END IF
   END DO
   END DO
   CALL MSPL_boxNpart_grad_t_deriv(s1,R,ndim,npart,Lbox,der_t4)
   PRINT *, "ALL-PARTICLES: ", SUM(der_t4), der_t4(s1%m/2,s1%Nknots/2,ndim/2,npart/2)
   PRINT *, "ONE-BY-ONE   : ", SUM(der_t5), der_t5(s1%m/2,s1%Nknots/2,ndim/2,npart/2)
   PRINT *, 

   PRINT *, "> > > Check the laplacian (MSPL_lapl_t_deriv)"
   ALLOCATE(der_t6(0:s1%m,0:s1%Nknots))
   ALLOCATE(lapl(0:s1%m,0:s1%Nknots))
   DO i2 = 1, 3, 1
      CALL RANDOM_NUMBER(x)
      j1=x*(s1%m+1)
      CALL RANDOM_NUMBER(x)
      j2=x*s1%Nknots
      CALL RANDOM_NUMBER(dist(1:ndim))
      dist(1:ndim)=s1%La+dist(1:ndim)*(s1%Lb-s1%La)*0.5d0
      dist(0)=DSQRT(DOT_PRODUCT(dist(1:ndim),dist(1:ndim)))
      !compute the function itself in der_t
      CALL MSPL_t_deriv(s1,dist(0),der_t)
      lapl=0.d0
      dv=0.0001d0
      DO i1 = 1, ndim, 1
         !forward, store in der_t2
         y=dist(i1)
         dist(i1)=dist(i1)+dv
         dist(0)=DSQRT(DOT_PRODUCT(dist(1:ndim),dist(1:ndim)))
         CALL MSPL_t_deriv(s1,dist(0),der_t2)
         dist(i1)=y
         dist(0)=DSQRT(DOT_PRODUCT(dist(1:ndim),dist(1:ndim)))
         !backward, store in der_t3(:,:,1)
         y=dist(i1)
         dist(i1)=dist(i1)-dv
         dist(0)=DSQRT(DOT_PRODUCT(dist(1:ndim),dist(1:ndim)))
         CALL MSPL_t_deriv(s1,dist(0),der_t3(:,:,1))
         dist(i1)=y
         dist(0)=DSQRT(DOT_PRODUCT(dist(1:ndim),dist(1:ndim)))
         !accumulate the laplacian
         lapl(:,:)=lapl(:,:)+(der_t2(:,:)-2.d0*der_t(:,:)+der_t3(:,:,1))/(dv*dv)
      END DO
      PRINT *, "Numerical gradient = ", SUM(lapl(:,:))
      CALL MSPL_lapl_t_deriv(s1,dist(0:ndim),ndim,der_t6)
      PRINT *, "Analitical laplacian = ", SUM(der_t6(:,:))
      PRINT *, 
   END DO

   PRINT *, "> > > Check the gradient_and_laplacian (MSPL_grad_and_lapl_t_deriv)"
   DO i2 = 1, 3, 1
      CALL RANDOM_NUMBER(x)
      j1=x*(s1%m+1)
      CALL RANDOM_NUMBER(x)
      j2=x*s1%Nknots
      CALL RANDOM_NUMBER(dist(1:ndim))
      dist(1:ndim)=s1%La+dist(1:ndim)*(s1%Lb-s1%La)*0.5d0
      dist(0)=DSQRT(DOT_PRODUCT(dist(1:ndim),dist(1:ndim)))
      CALL MSPL_grad_t_deriv(s1,dist(0:ndim),ndim,grad)
      PRINT *, "Gradient alone = ", SUM(grad(:,:,:))
      CALL MSPL_lapl_t_deriv(s1,dist(0:ndim),ndim,der_t6)
      PRINT *, "Laplacian alone = ", SUM(der_t6(:,:))
      CALL MSPL_grad_and_lapl_t_deriv(s1,dist(0:ndim),ndim,grad,der_t6)
      PRINT *, "Gradient together with laplacian = ", SUM(grad)
      PRINT *, "Laplacian together with gradient = ", SUM(der_t6)
      PRINT *, 
   END DO

   PRINT *, "> > > Check the all-particles gradient_and_laplacian (MSPL_boxNpart_grad_and_lapl_t_deriv)"
   ALLOCATE(der_t7(0:s1%m,0:s1%Nknots,1:npart))
   CALL RANDOM_NUMBER(R)
   R=s1%La+R*(s1%Lb-s1%La)
   Lbox=s1%Lb-s1%La
   CALL MSPL_boxNpart_grad_and_lapl_t_deriv(s1,R,ndim,npart,Lbox,der_t5,der_t7)
   PRINT *, "ALL-PARTICLES GRADIENT: ", SUM(der_t5), der_t5(s1%m/2,s1%Nknots/2,ndim/2,npart/2)
   PRINT *, "ALL-PARTICLES LAPLACIAN: ", SUM(der_t7), der_t7(s1%m/2,s1%Nknots/2,npart/2)
   der_t7=0.d0
   der_t5=0.d0
   DO i2 = 1, npart, 1
   DO i1 = 1, npart, 1
      IF (i1 /= i2) THEN
         dist(1:ndim)=R(1:ndim,i2)-R(1:ndim,i1)
         dist(1:ndim)=dist(1:ndim)-Lbox(1:ndim)*DNINT(dist(1:ndim)/Lbox(1:ndim))
         dist(0)=DSQRT(DOT_PRODUCT(dist(1:ndim),dist(1:ndim)))
         CALL MSPL_grad_and_lapl_t_deriv(s1,dist(0:ndim),ndim,der_t4(:,:,:,1),lapl)
         der_t5(:,:,:,i2)=der_t5(:,:,:,i2)+der_t4(:,:,:,1)
         der_t7(:,:,i2)=der_t7(:,:,i2)+lapl
      END IF
   END DO
   END DO
   PRINT *, "ONE-BY-ONE GRADIENT: ", SUM(der_t5), der_t5(s1%m/2,s1%Nknots/2,ndim/2,npart/2)
   PRINT *, "ONE-BY-ONE LAPLACIAN: ", SUM(der_t7), der_t7(s1%m/2,s1%Nknots/2,npart/2)
   PRINT *, 
   
   DEALLOCATE(der_t7)
   DEALLOCATE(lapl)
   DEALLOCATE(der_t6)
   DEALLOCATE(der_t5)
   DEALLOCATE(der_t4)
   DEALLOCATE(grad)
   DEALLOCATE(der_t3)
   DEALLOCATE(der_t2)
	DEALLOCATE(der_t)
	
	CALL MSPL_deallocate(s1)
	CALL MSPL_deallocate(s2)

END PROGRAM prova
