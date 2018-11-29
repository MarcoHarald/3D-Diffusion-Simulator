PROGRAM pdes
IMPLICIT NONE

INTEGER,PARAMETER :: n=300
INTEGER :: i, j, kkk
REAL :: x,y,tolerance,accuracy,temp(-n:n,-n:n),dx_temp(-n:n,-n:n)
REAL :: c_dfsn,c_r,d_r,d_t, temp0, temp1
!SET VALUES

	

!diffusion constant (cm2 s-1)
!radius of initial (cm)
!step size (cm)
!step size (s)

	c_dfsn = 1.1
	c_r = 1.0
	d_r = 0.1
	d_t = 1.0e-4

	!CONNECT FILES FOR INPUT/OUTPUT
	OPEN(UNIT=12, FILE="pde.dat")
        

!***********************
WRITE(6,*) 'CHECK'
!***********************

!SET ACCURACY & TOLERANCE FOR BOUNDARY LIMITS
tolerance = 10.0**(-5.0)
accuracy = 10.0**(-3.0)

!SET INITIAL GUESS AT BOUNDARY CONDITIONS (BC)
temp1 = 273.0
temp0 = 293.0


!***********************
WRITE(6,*) 'CHECK'
!***********************


!CHECK IF COORDINATE IS INSIDE DNA CYLINDER & SET IT TO BC
DO i=-n,n
	DO j=-n,n
		x=REAL(i)*d_r
		y=REAL(j)*d_r
		IF(x**2.0+y**2.0<c_r**2.0+tolerance)temp(i,j)=temp0
	END DO
END DO


!***********************
WRITE(6,*) 'CHECK'
!***********************


!LOOP OVER ALL COORDINATES
do kkk=1,100
	do i=-n+1,n-1
		do j=-n+1,n-1

			!ASSIGN ACTUAL COORDINATE VALUES TO VARS
			x=real(i)*d_r
			y=real(j)*d_r
			!IF INSIDE CYLINDER, NO ACTION
			if(x**2.0+y**2.0<c_r**2.0+tolerance)then

			else
				!ESTIMATION ITERATOR
					dx_temp(i) = (temp(i+1)-2.0*temp(i)+temp(i-1))/d_r**2.0				
				!FOR CHANGE BIGGER THAN accuracy MARGIN, THEN INITIATE ANOTHER ITERATION

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
			x=real(i)*d_r
			y=real(j)*d_r
			!WRITE TO FILE
			WRITE(12,*) x, y ,dx_temp(i,j)

		enddo
	enddo


!***********************
WRITE(6,*) 'CHECK'
!***********************



END PROGRAM pdes

