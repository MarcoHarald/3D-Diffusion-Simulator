!evaluating the effect of temperature hitting the boundaries


PROGRAM pdes
IMPLICIT NONE

INTEGER,PARAMETER :: n=30, m=30
INTEGER :: i, j, k,l, steps
REAL :: x(-n:n),y(-m:m),tolerance,temp(-n:n,-m:m),dx_temp(-n:n),dy_temp(-m:m),dt_temp(-n:n, -m:m)    ! array(-n:n,0:100)
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


	c_dfsn = 1.1
	c_r = 1.0
	d_r = 0.1
	d_t = 1.0e-0

	!CONNECT FILES FOR INPUT/OUTPUT
	OPEN(UNIT=12, FILE="dif.dat")
 	OPEN(UNIT=13, FILE="dmn.dat")
	OPEN(UNIT=14, FILE="dif2.dat")
!***********************
WRITE(6,*) 'CHECK   1'
!***********************

!SET ACCURACY & TOLERANCE FOR BOUNDARY LIMITS
tolerance = 10.0**(-5.0)

!SET INITIAL GUESS AT BOUNDARY CONDITIONS (BC)
temp1 = 0.0
temp0 = 20.0

!SET SIMULATION LENGTH (s)
t_siml = 100.0

!DEFINE EXACT VALUE OF COORDINATES ON TWO AXES
DO i=-n,n
x(i)=REAL(i)*d_r
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
	 
ALLOCATE(array(-x_position:x_position,-y_position:y_position,t_frame)) 
!DEALLOCATE(array)


!***********************
WRITE(6,*) 'CHECK   2'
!***********************


!CHECK IF COORDINATE IS INSIDE INITIAL CONDITION & SET IT TO BC
!INITIALISE FIRST 'ARRAY' VALUES (INITIAL CONDITIONS)
DO j=-m,m	
	DO i=-n,n
		IF(x(i)**2.0<c_r**2.0+tolerance .AND. y(j)**2.0<c_r**2.0+tolerance) THEN
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
!OUTPUT:end values (loop cycles over x-position & then y-position)
DO l=-m,m	
	DO i=-n,n	
		WRITE(12,*) x(i), y(l),array(i,j,1)
	END DO 

	WRITE(12,*)
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
			dx_temp(j) = ((array(j+1,l,1))-(2.0*array(j,l,1))+(array(j-1,l,1)))/(d_r**2.0)
			dy_temp(l) = ((array(j,l+1,1))-(2.0*array(j,l,1))+(array(j,l-1,1)))/(d_r**2.0)	
		END DO
	END DO

!RATE OF CHANGE OF TEMP WITH RESPECT TO TIME & ITERATE TO ACCOUNT FOR CHANGE IN TEMP
	DO l=-m+1, m-1	
		DO i=-n+1,n-1
			dt_temp(i,l) = c_dfsn*dx_temp(i)+c_dfsn*dy_temp(l)
			array(i,l,k) = array(i,l,k-1) + (d_t*dt_temp(i,l))
		END DO
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


!OUTPUT:end values (loop cycles over x-position & then y-position)
DO l=-m,m	
	DO i=-n,n	
		WRITE(14,*) x(i), y(l),array(i,l,steps-1)
	END DO 

	WRITE(14,*)
END DO

!OUTPUT:temperature distribution over time (loop 'samples' the data every 1000th value to save memory)
	DO k=1,steps,1000
		DO  i=-n,n
			WRITE(13,*) x(i),k, array(i,0,k) 		
		END DO
!leave gap
		WRITE(13,*)
	END DO

!***********************
WRITE(6,*) 'CHECK   6'
!***********************



END PROGRAM pdes

