PROGRAM pdes
IMPLICIT NONE

INTEGER,PARAMETER :: n=80
INTEGER :: i, j, flag, kkk
REAL :: x,y, phi_s,tolerance, kappa_2, c, e, eps, k, temp, h, r_dna, accuracy
REAL :: phi(-n:n,-n:n), oldval

!SET VALUES
!speed of light 			c
!e constant						e
!epsilon							eps
!boltzmann constant		k
!temperature					t
!radius 							r_dna
!step size						h

	c = 0.3
	e = 1.6E-19
	eps = 80.0*8.85e-12
	k = 1.38064852*(10.0**(-23.0))
	temp = 310.0
	!r_dna = 1.0
	!h = r_dna/real(n)
	h = 0.05 


	!CONNECT FILES FOR INPUT/OUTPUT
	OPEN(UNIT=12, FILE="pde.dat")
        

!***********************
WRITE(6,*) 'CHECK'
!***********************

!SET ACCURACY & TOLERANCE FOR BOUNDARY LIMITS
tolerance = 10.0**(-5.0)
accuracy = 10.0**(-3.0)

!SET INITIAL GUESS AT BOUNDARY CONDITIONS (BC)
phi = 0.0
phi_s=0.1

!SET DEBYE LENGTH
kappa_2 = (2.0*c*e**2.0)/(eps*k*temp)

!***********************
WRITE(6,*) 'CHECK'
!***********************


!CHECK IF COORDINATE IS INSIDE DNA CYLINDER & SET IT TO BC
DO i=-n,n
	DO j=-n,n
		x=REAL(i)*h
		y=REAL(j)*h
		IF(x**2.0+y**2.0<r_dna**2.0+tolerance)phi(i,j)=phi_s
	END DO
END DO


!***********************
WRITE(6,*) 'CHECK'
!***********************


!LOOP OVER ALL COORDINATES
do kkk=1,100
	do i=-n+1,n-1
		do j=-n+1,n-1

		!	oldval=phi(i,j)
			!ASSIGN ACTUAL COORDINATE VALUES TO VARS
			x=real(i)*h
			y=real(j)*h
			!IF INSIDE CYLINDER, NO ACTION
			if(x**2.0+y**2.0<r_dna**2.0+tolerance)then

			else
				!ESTIMATION ITERATOR
				phi(i,j)= (phi(i+1,j)+phi(i-1,j)+phi(i,j+1)+phi(i,j-1))/(4.0+h**2.0*kappa_2)
				!FOR CHANGE BIGGER THAN accuracy MARGIN, THEN INITIATE ANOTHER ITERATION
			!	if(abs(oldval-phi(i,j)) < accuracy) flag=1	
			endif
		enddo
	enddo

END DO


!***********************
WRITE(6,*) 'CHECK'
!***********************




	do i=-n+1,n-1
		do j=-n+1,n-1

			!ASSIGN ACTUAL COORDINATE VALUES TO VARS
			x=real(i)*h
			y=real(j)*h
			!WRITE TO FILE
			WRITE(12,*) x, y ,phi(i,j)

		enddo
	enddo


!***********************
WRITE(6,*) 'CHECK'
!***********************



END PROGRAM pdes

