!evaluating the effect of temperature hitting the boundaries


PROGRAM pdes
IMPLICIT NONE

INTEGER,PARAMETER :: n=300
INTEGER :: i, j, k, steps
REAL :: x(-n:n),tolerance,temp(-n:n),dx_temp(-n:n),dt_temp(-n:n)    ! array(-n:n,0:100)
REAL :: c_dfsn,c_r,d_r,d_t, temp0, temp1,t_siml	


!2-dimensional array of INTEGERS with size to be allocated
REAL, DIMENSION(:,:), ALLOCATABLE :: array
!allocatable integers to define dimensions
INTEGER :: x_position, t_frame
	

!----------------REFERENCE OF VARIABLES-----------------

!array		motion array over time (save slices to represent over time)
!x_position	point along the defined space (on which a temp is assigned)
!t_frame		snapshot of instance in time


!c_dfsn		diffusion constant (cm2 s-1)
!c_r		radius of initial (cm)
!d_r		space step size (cm)
!d_t		time step size (s)

!	***ACTION***	change variable names for demonstrators!
!			full width half maximum

!----------------END OF REFERENCE-----------------


	c_dfsn = 1.1
	c_r = 1.0
	d_r = 0.1
	d_t = 1.0e-4

	!CONNECT FILES FOR INPUT/OUTPUT
	OPEN(UNIT=12, FILE="dif.dat")
 	OPEN(UNIT=13, FILE="dmn.dat")

!***********************
WRITE(6,*) 'CHECK   1'
!***********************

!SET ACCURACY & TOLERANCE FOR BOUNDARY LIMITS
tolerance = 10.0**(-5.0)

!SET INITIAL GUESS AT BOUNDARY CONDITIONS (BC)
temp1 = 0.0
temp0 = 20.0

!SET SIMULATION LENGTH (s)
t_siml = 5.0


!INITIATE STEPS VARIABLE
steps = int(real(t_siml)/d_t)
WRITE(6,*) '# OF STEPS: ', steps

!DETERMINE ALLOCATABLE ARRAY DIMENSIONS
	x_position = 300
	t_frame = steps
	 
ALLOCATE(array(-x_position:x_position,t_frame)) 
!DEALLOCATE(array)


!***********************
WRITE(6,*) 'CHECK   2'
!***********************


!CHECK IF COORDINATE IS INSIDE INITIAL CONDITION & SET IT TO BC
!INITIALISE FIRST 'ARRAY' VALUES (INITIAL CONDITIONS)
DO i=-n,n
	x(i)=REAL(i)*d_r
	IF(x(i)**2.0<c_r**2.0+tolerance) THEN
			temp(i)=temp0
			array(i, 1)=temp(i)
	ELSEIF(x(i)>-29.0 .AND. x(i)<-27.0+tolerance) THEN
			temp(i)=temp0
			array(i, 1)=temp(i)

		ELSE
			temp(i)=temp1
			array(i, 1)=temp(i)
	END IF 
END DO


!***********************
WRITE(6,*) 'CHECK  3'
!***********************


!LOOP OVER INSTANCES OF TIME (FRAMES)
DO k=1, steps

!SECOND ORDER DIFF IN SPACE
	DO j=(-n+1),n-1
		dx_temp(j) = ((temp(j+1))-(2.0*temp(j))+(temp(j-1)))/(d_r**2.0)									
	END DO

!RATE OF CHANGE OF TEMP WITH RESPECT TO TIME & ITERATE TO ACCOUNT FOR CHANGE IN TEMP
	DO i=-n+1,n-1
		dt_temp(i) = c_dfsn*dx_temp(i)
		temp(i) = temp(i) + (d_t*dt_temp(i))
		array(i,k) = temp(i) 
	END DO

END DO

!***********************
WRITE(6,*) 'CHECK   4A '
!***********************

!END DO
!FOR CHANGE BIGGER THAN accuracy MARGIN, THEN INITIATE ANOTHER ITERATION

!***********************
WRITE(6,*) 'CHECK   5'
!***********************



!OUTPUT:end values (loop cycles over x-position)
	DO i=-n,n	
		WRITE(12,*) x(i),temp(i)
		WRITE(12,*)
	END DO 

!OUTPUT:temperature distribution over time (loop 'samples' the data every 1000th value to save memory)
	DO k=1,steps,1000
		DO  i=-n,n
			WRITE(13,*) x(i), k, array(i,k) 		
		END DO

		WRITE(13,*)
	END DO

!***********************
WRITE(6,*) 'CHECK   6'
!***********************



END PROGRAM pdes

