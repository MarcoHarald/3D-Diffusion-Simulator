!evaluating the effect of temperature hitting the boundaries


PROGRAM pdes
IMPLICIT NONE

INTEGER,PARAMETER :: n=10, m=10
INTEGER :: i, j, k,l, steps
REAL :: x(-n:n),y(-m:m),tolerance,temp(-n:n,-m:m),dx_temp(-n:n,-m:m),dy_temp(-n:n,-m:m),dt_temp(-n:n,-m:m)    ! array(-n:n,0:100)
REAL :: c_dfsn,c_r,d_r,d_t, temp0, temp1,t_siml	

CHARACTER :: checkBreak

!2-dimensional array of INTEGERS with size to be allocated
REAL, DIMENSION(:,:,:), ALLOCATABLE :: array
!allocatable integers to define dimensions
INTEGER :: x_position,y_position, t_frame
	

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


	c_dfsn = 1.5
	c_r = 3.0
	d_r = 1.0
	d_t = 1.0e-1

	!CONNECT FILES FOR INPUT/OUTPUT
	OPEN(UNIT=12, FILE="prelim.dat")
 	OPEN(UNIT=13, FILE="sampling.dat")
	OPEN(UNIT=14, FILE="endVal.dat")
	OPEN(UNIT=15, FILE="diffs.dat")
!	OPEN(UNIT=16, FILE="loopPrint.dat")
	OPEN(UNIT=22, FILE="endTest.dat")
!***********************
WRITE(6,*) 'CHECK   1'
!***********************

!SET ACCURACY & TOLERANCE FOR BOUNDARY LIMITS
tolerance = 10.0**(-5.0)

!SET INITIAL GUESS AT BOUNDARY CONDITIONS (BC)
temp1 = 10.0
temp0 = 20.0

!SET SIMULATION LENGTH (s)
t_siml = 10.0

!DEFINE EXACT VALUE OF COORDINATES ON TWO AXES
DO i=-n,n
x(i)=REAL(i)*d_r
WRITE(6,*) x(i)
END DO

DO i=-m,m
y(i)=REAL(i)*d_r
END DO

!INITIATE STEPS VARIABLE
steps = int(real(t_siml)/d_t)
WRITE(6,*) '# OF STEPS: ', steps

!DETERMINE ALLOCATABLE ARRAY DIMENSIONS
	x_position = n
	y_position = m
	t_frame = steps

!----------------------------------------------
WRITE(6,*)'X-AXIS WIDTH: ', n
WRITE(6,*)'Y-AXIS WIDTH: ', m
!----------------------------------------------

ALLOCATE(array(-x_position:x_position,-y_position:y_position,1:t_frame)) 
!DEALLOCATE(array)


!***********************
WRITE(6,*) 'CHECK   2'
!***********************


!CHECK IF COORDINATE IS INSIDE INITIAL CONDITION & SET IT TO BC
!INITIALISE FIRST 'ARRAY' VALUES (INITIAL CONDITIONS)
DO j=-m,m	
	DO i=-n,n
		IF(abs(x(i))<c_r+tolerance .AND. abs(y(j))<c_r+tolerance) THEN
				array(i, j, 1)=temp0
		!ELSEIF(x(i)>-29.0 .AND. x(i)<-27.0+tolerance) THEN
		!		temp(i)=temp0
		!		array(i, j, 1)=temp(i)

			ELSE
				array(i, j, 1)=temp1
		END IF 
	END DO
END DO

!***********************
WRITE(6,*) 'CHECK  3'
!***********************
!OUTPUT:initial values (loop cycles over x-position & then y-position)
DO l=-m,m	
	DO i=-n,n	
		WRITE(12,*) x(i), y(l),array(i,l,1)
	END DO 

!	WRITE(12,*)
END DO
!-----------------------
!WRITE(6,*) k,'b'
!WRITE(6,*) '(d_r**2.0)',(d_r**2.0)
!WRITE(6,*) '(temp(j+1,l)', (temp(j+1,l))
!WRITE(6,*) '(2.0*temp(j,l))', (2.0*temp(j,l))
!WRITE(6,*) '(temp(j-1,l))', (temp(j-1,l))
!WRITE(6,*) 
!WRITE(6,*)
!------------------------
!------------------------
!WRITE(6,*) k,'c'
!------------------------
!------------------------
!WRITE(6,*) 'Continue?   [INSERT CHARACTER]'
!READ(5,*)checkBreak	
!------------------------


!LOOP OVER INSTANCES OF TIME (FRAMES)
DO k=2, steps
!GRADIENT OF TEMP CHANGE WITH RESPECT TO TIME [ x ]
	DO l=(-m+1),m-1
		DO j=(-n+1),n-1
			dx_temp(j,l) = ((array(j+1,l,k))-(2.0*array(j,l,k))+(array(j-1,l,k)))/(d_r*d_r)
			dy_temp(j,l) = ((array(j,l+1,k))-(2.0*array(j,l,k))+(array(j,l-1,k)))/(d_r*d_r)	
			
			!-----------------------
			!WRITE(16,*) k,array(j+1,l,1),(2.0*array(j,l,1)),array(j-1,l,1)
			!-----------------------
					
			!-----------------------
				IF(dx_temp(j,l)>tolerance .OR. dy_temp(j,l)>tolerance) THEN
					WRITE(15,*) k, dx_temp,dy_temp			
				ENDIF
			!-----------------------
		
		END DO
	END DO

	!-----------------------
	!WRITE(16,*)
	!-----------------------
	!-----------------------
	!WRITE(15,*)
	!-----------------------


!RATE OF CHANGE OF TEMP WITH RESPECT TO TIME & ITERATE TO ACCOUNT FOR CHANGE IN TEMP
	DO l=-m+1, m-1	
		DO i=-n+1,n-1
			dt_temp(i,l) = c_dfsn*dx_temp(i,l)+c_dfsn*dy_temp(i,l)
			array(i,l,k) = array(i,l,k-1) + (d_t*dt_temp(i,l))
		END DO
	END DO
END DO

!***********************
WRITE(6,*) 'CHECK   4'
!***********************

!END DO
!FOR CHANGE BIGGER THAN accuracy MARGIN, THEN INITIATE ANOTHER ITERATION

!***********************
WRITE(6,*) 'CHECK   5'
!***********************


!OUTPUT:end values (loop cycles over x-position & then y-position)
DO l=1-m,m-1	
	DO i=1-n,n-1	
		WRITE(14,*) x(i), y(l),array(i,l,steps-1)
	END DO 

	WRITE(14,*)
END DO

!OUTPUT:temperature distribution over time (loop 'samples' the data every 1000th value to save memory)
	DO k=1,steps,1000
		DO  i=1-n,n-1
			WRITE(13,*) x(i),k, array(i,0,k) 		
		END DO
!leave gap
		WRITE(13,*)
	END DO

!***********************
WRITE(6,*) 'CHECK   6'
!***********************


!--------------------------------------------------
!OUTPUT:initial values (loop cycles over x-position & then y-position)
k=1
DO l=0, 0	
	DO i=1-n,n-1	
		WRITE(22,*) x(i), y(l),array(i,l,k)
	END DO 

	WRITE(22,*)
END DO
!--------------------------------------------------


END PROGRAM pdes

