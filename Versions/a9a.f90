!evaluating the effect of temperature hitting the boundaries


PROGRAM pdes
IMPLICIT NONE

INTEGER,PARAMETER :: x0=10, y0=10, z0=10
INTEGER :: x,y,z,i,k, steps
REAL :: xAxis(-x0:x0),yAxis(-y0:y0),zAxis(-z0:z0)
REAL :: dx_temp(-x0:x0,-y0:y0,-z0:z0),dy_temp(-x0:x0,-y0:y0,-z0:z0),dz_temp(-x0:x0,-y0:y0,-z0:z0)
REAL :: tolerance, dt_temp(-x0:x0,-y0:y0,-z0:z0) !temp(-x0:x0,-y0:y0,-z0:z0)
REAL :: c_dfsn,c_r,d_r,d_t, temp0, temp1,t_siml	

!4-dimensional array of REAL with size to be allocated
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: array
!allocatable integers to define dimensions
INTEGER :: x_position,y_position,z_position, t_frame, stepSize
	

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
	c_r = 2.0
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

!DEFINE EXACT VALUE OF COORDINATES ON THREE AXES
DO i=-x0,x0
xAxis(i)=REAL(i)*d_r
END DO

DO i=-y0,y0
yAxis(i)=REAL(i)*d_r
END DO

DO i=-z0,z0
zAxis(i)=REAL(i)*d_r
END DO

!INITIATE STEPS VARIABLE
steps = int(real(t_siml)/d_t)
WRITE(6,*) '# OF STEPS: ', steps

!DETERMINE ALLOCATABLE ARRAY DIMENSIONS
	x_position = x0
	y_position = y0
	z_position = z0
	t_frame = steps

!----------------------------------------------
WRITE(6,*)'X-AXIS WIDTH: ', x0
WRITE(6,*)'Y-AXIS WIDTH: ', y0
WRITE(6,*)'Z-AXIS WIDTH: ', z0
!----------------------------------------------

ALLOCATE(array(-x_position:x_position,-y_position:y_position,-z_position:z_position,1:t_frame)) 
!DEALLOCATE(array)


!***********************
WRITE(6,*) 'CHECK   2'
!***********************


!CHECK IF COORDINATE IS INSIDE INITIAL CONDITION & SET IT TO BC
!INITIALISE FIRST 'ARRAY' VALUES (INITIAL CONDITIONS)
DO z=-z0,z0
	DO y=-y0,y0	
		DO x=-x0,x0
			IF(abs(xAxis(x))<c_r+tolerance .AND. abs(yAxis(y))<c_r+tolerance .AND. abs(zAxis(z))<c_r+tolerance) THEN
					array(x,y,z,1)=temp0
				ELSE
					array(x,y,z,1)=temp1
			END IF 
		END DO
	END DO
END DO

!***********************
WRITE(6,*) 'CHECK  3'
!***********************


!OUTPUT:initial values (loop cycles over x-position & then y-position)
DO z=-z0,z0
	DO y=-y0,y0	
		DO x=-x0,x0	
		WRITE(12,*) xAxis(x),yAxis(y),zAxis(z),array(x,y,z,1)
		END DO 
	END DO
	WRITE(12,*)
END DO


!***********************
WRITE(6,*) 'CHECK  3A'
!***********************

!LOOP OVER INSTANCES OF TIME (FRAMES)
DO k=1, steps-1
!GRADIENT OF TEMP CHANGE WITH RESPECT TO TIME [ x ]
	DO z=1-z0,z0-1
		DO y=1-y0,y0-1	
			DO x=1-x0,x0-1
				dx_temp(x,y,z) = ((array(x+1,y,z,k))-(2.0*array(x,y,z,k))+(array(x-1,y,z,k)))/(d_r*d_r)
				dy_temp(x,y,z) = ((array(x,y+1,z,k))-(2.0*array(x,y,z,k))+(array(x,y-1,z,k)))/(d_r*d_r)
				dZ_temp(x,y,z) = ((array(x,y,z+1,k))-(2.0*array(x,y,z,k))+(array(x,y,z-1,k)))/(d_r*d_r)			
			END DO
		END DO
	END DO

!RATE OF CHANGE OF TEMP WITH RESPECT TO TIME & ITERATE TO ACCOUNT FOR CHANGE IN TEMP
	DO z=1-z0,z0-1
		DO y=1-y0,y0-1	
			DO x=1-x0,x0-1
				dt_temp(x,y,z) = c_dfsn*dx_temp(x,y,z)+c_dfsn*dy_temp(x,y,z)+c_dfsn*dz_temp(x,y,z)
				array(x,y,z,k+1) = array(x,y,z,k) + (d_t*dt_temp(x,y,z))
			END DO
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


!OUTPUT:end values for all 3D space (loop cycles over x-position & then y-position)
	DO z=1-z0,z0-1
		DO y=1-y0,y0-1	
			DO x=1-x0,x0-1
				WRITE(14,*) xAxis(x), yAxis(y),zAxis(z),array(x,y,z,steps-1)
			END DO 
		END DO
		WRITE(14,*)
	END DO

!OUTPUT:temperature distribution for a slice of the 3D space, over time (loop 'samples' the data every 1000th value to save memory)
	DO k=1,steps,1000
		DO  x=1-x0,x0-1
			WRITE(13,*) xAxis(x),yAxis(0),zAxis(0),k,array(x,0,0,k) 		
		END DO
		WRITE(13,*)
	END DO

!***********************
WRITE(6,*) 'CHECK   6'
!***********************


!--------------------------------------------------
!OUTPUT:slice of final values (loop cycles slices of the z axis)
k=steps
stepSize = 5

	DO i=0,2*z0-stepSize, stepSize
		z=-z0+1+i

		OPEN(UNIT=23+i/stepSize, FILE=trim("endVal")//trim(char(i))//".dat")

		DO y=1-y0,y0-1	
			DO x=1-x0,x0-1
				WRITE(22+i,*) xAxis(x),yAxis(y),array(x,y,z,k)
			END DO 
		END DO
		WRITE(22,*)
	END DO
!--------------------------------------------------


END PROGRAM pdes

