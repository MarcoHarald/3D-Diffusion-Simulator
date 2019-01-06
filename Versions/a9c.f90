PROGRAM tempDiffusion
IMPLICIT NONE

INTEGER, PARAMETER :: dP = SELECTED_REAL_KIND(15)			!Set double precision
INTEGER,PARAMETER :: x0=10, y0=10, z0=10					!Set simulator space dimensions
INTEGER :: x,y,z,i,k, steps, flag, boundInt, lowVal, highVal
REAL(dP) :: xAxis(-x0:x0),yAxis(-y0:y0),zAxis(-z0:z0)
REAL(dP) :: dx_temp(-x0:x0,-y0:y0,-z0:z0),dy_temp(-x0:x0,-y0:y0,-z0:z0),dz_temp(-x0:x0,-y0:y0,-z0:z0)
REAL(dP) :: tolerance, dt_temp(-x0:x0,-y0:y0,-z0:z0)
REAL(dP) :: constDfsn,initRad,stepDist,stepTime, temp0, temp1,t_siml, bound, initRad2, tempPeak

!4-dimensional array of REAL with size to be allocated
REAL(dP), DIMENSION(:,:,:,:), ALLOCATABLE :: array
REAL(dP), DIMENSION(:), ALLOCATABLE :: gaussian, fwhm

!allocatable integers to define dimensions
INTEGER :: x_position,y_position,z_position, t_frame, stepSize
	

!----------------REFERENCE OF VARIABLES-----------------

!array		motion array over time (save slices to represent over time)
!x_position	point along the defined space (on which a temp is assigned)
!t_frame		snapshot of instance in time


!constDfsn		constant of diffusion  (cm2 s-1)
!initRad		initial radius (cm)
!stepDist		space step size (cm)
!stepTime		time step size (s)

!	***ACTION***	change variable names for demonstrators
!			full width half maximum
!			evaluate the effect of temperature hitting the boundaries


!----------------END OF REFERENCE-----------------


	constDfsn = 1.1
	initRad = 0.2
	initRad2 = 4.0
	stepDist = 1.0
	stepTime = 1.0e-1

	!CONNECT FILES FOR INPUT/OUTPUT
	OPEN(UNIT=12, FILE="initAll.dat")
	OPEN(UNIT=13, FILE="endAll.dat")
	OPEN(UNIT=14, FILE="initSlice.dat")
    OPEN(UNIT=15, FILE="endSlice.dat")
    OPEN(UNIT=16, FILE="fwhm.dat")
	OPEN(UNIT=30, FILE="val0.dat")
	OPEN(UNIT=31, FILE="val1.dat")
	OPEN(UNIT=32, FILE="val2.dat")
	OPEN(UNIT=33, FILE="val3.dat")

!***********************
WRITE(6,*) 'CHECK   1'
!***********************

!SET ACCURACY & TOLERANCE FOR BOUNDARY LIMITS
tolerance = 10.0**(-5.0)

!SET INITIAL GUESS AT BOUNDARY CONDITIONS (BC)
temp1 = 0.0
temp0 = 20.0

!SET SIMULATION LENGTH (s)
t_siml = 10.0

!DEFINE EXACT VALUE OF COORDINATES ON THREE AXES
DO i=-x0,x0
	xAxis(i)=REAL(i)*stepDist
END DO

DO i=-y0,y0
	yAxis(i)=REAL(i)*stepDist
END DO

DO i=-z0,z0
	zAxis(i)=REAL(i)*stepDist
END DO

!INITIATE STEPS VARIABLE
steps = int(real(t_siml)/stepTime)
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

ALLOCATE(array(-x_position:x_position,-y_position:y_position,-z_position:z_position,0:t_frame)) 
!DEALLOCATE(array)

!***********************
WRITE(6,*) 'CHECK   2'
!***********************

!CHECK IF COORDINATE IS INSIDE INITIAL CONDITION & SET IT TO BC
!INITIALISE FIRST 'ARRAY' VALUES (INITIAL CONDITIONS)
DO z=-z0,z0
	DO y=-y0,y0	
		DO x=-x0,x0
			IF((xAxis(x)**2+yAxis(y)**2+zAxis(z)**2)<initRad**2) THEN
				array(x,y,z,0)=temp0
			ELSE
				array(x,y,z,0)=temp1
			ENDIF 
		END DO
	END DO
END DO

!***********************
WRITE(6,*) 'CHECK  3'
!***********************

!LOOP OVER INSTANCES OF TIME (FRAMES) (set to steps-1] to avoid bug)
DO k=0, steps-1
!GRADIENT OF TEMP CHANGE WITH RESPECT TO TIME [ x ]
	DO z=1-z0,z0-1
		DO y=1-y0,y0-1	
			DO x=1-x0,x0-1
				dx_temp(x,y,z) = (array(x+1,y,z,k)-(2.0*array(x,y,z,k))+array(x-1,y,z,k))/(stepDist**2)
				dy_temp(x,y,z) = (array(x,y+1,z,k)-(2.0*array(x,y,z,k))+array(x,y-1,z,k))/(stepDist**2)
				dz_temp(x,y,z) = (array(x,y,z+1,k)-(2.0*array(x,y,z,k))+array(x,y,z-1,k))/(stepDist**2)				
			END DO
		END DO
	END DO

!RATE OF CHANGE OF TEMP WITH RESPECT TO TIME & ITERATE TO ACCOUNT FOR CHANGE IN TEMP
	DO z=1-z0,z0-1
		DO y=1-y0,y0-1	
			DO x=1-x0,x0-1
                dt_temp(x,y,z) = constDfsn*(dx_temp(x,y,z)+dy_temp(x,y,z)+dz_temp(x,y,z))
                
				!to eliminate any binary number errors, only values above a threshold are accepted (above the double precision number precision of 15 significant figures)
				IF(abs(dt_temp(x,y,z))>1.0E-8) THEN
					array(x,y,z,k+1) = array(x,y,z,k) + (stepTime*dt_temp(x,y,z))
				ELSE
					array(x,y,z,k+1) = array(x,y,z,k)
					dt_temp(x,y,z) = 0
				ENDIF
			END DO
		END DO
	END DO

!ADJUST TEMPERATURE AT EDGES OF SIMULATION-SPACE TO EMULATE AN INFINITELY LARGE SPACE
	!done for all six planes of the parallelepipid that encompasses the simulation space
	!simple estimate of the temperature change trend of two adjacent inner facing blocks
	x=-x0
	DO y=1-y0,y0-1	
		DO x=1-x0,x0-1
			IF (abs(dt_temp(x+1,y,z))>= 0+tolerance) THEN
				dt_temp(x,y,z) = dt_temp(x+1,y,z)+(dt_temp(x+1,y,z)- dt_temp(x+2,y,z)) 
				array(x,y,z,k+1) = array(x,y,z,k) + (stepTime*dt_temp(x,y,z))
			ENDIF
		END DO
	END DO

	x=x0
	DO y=1-y0,y0-1	
		DO x=1-x0,x0-1
			IF (abs(dt_temp(x-1,y,z))>= 0+tolerance) THEN
				dt_temp(x,y,z) = dt_temp(x-1,y,z)+(dt_temp(x-1,y,z)- dt_temp(x-2,y,z)) 
				array(x,y,z,k+1) = array(x,y,z,k) + (stepTime*dt_temp(x,y,z))
			ENDIF	
		END DO
	END DO	

	y=-y0
		DO z=1-z0,z0-1	
			DO x=1-x0,x0-1
				IF (abs(dt_temp(x,y+1,z))>= 0+tolerance) THEN
					dt_temp(x,y,z) = dt_temp(x,y+1,z)+(dt_temp(x,y+1,z)- dt_temp(x,y+2,z)) 
					array(x,y,z,k+1) = array(x,y,z,k) + (stepTime*dt_temp(x,y,z))
				ENDIF
			END DO
		END DO

	y=y0
	DO z=1-z0,z0-1	
		DO x=1-x0,x0-1
			IF (abs(dt_temp(x,y-1,z))>= 0+tolerance) THEN
				dt_temp(x,y,z) = dt_temp(x,y-1,z)+(dt_temp(x,y-1,z)- dt_temp(x,y-2,z)) 
				array(x,y,z,k+1) = array(x,y,z,k) + (stepTime*dt_temp(x,y,z))
			END IF
		END DO
	END DO	

	z=-z0
		DO y=1-y0,y0-1	
			DO x=1-x0,x0-1
				IF (abs(dt_temp(x,y,z+1))>= 0+tolerance) THEN
					dt_temp(x,y,z) = dt_temp(x,y,z+1)+(dt_temp(x,y,z+1)- dt_temp(x,y,z+2)) 
					array(x,y,z,k+1) = array(x,y,z,k) + (stepTime*dt_temp(x,y,z))
				ENDIF		
			END DO
		END DO

    
        z=z0
	DO y=1-y0,y0-1	
		DO x=1-x0,x0-1
			IF (abs(dt_temp(x,y,z-1))>= 0+tolerance) THEN
				dt_temp(x,y,z) = dt_temp(x,y,z-1)+(dt_temp(x,y,z-1)- dt_temp(x,y,z-2)) 
				array(x,y,z,k+1) = array(x,y,z,k) + (stepTime*dt_temp(x,y,z))
			ENDIF
		END DO
	END DO	
END DO

!***********************
WRITE(6,*) 'CHECK   4'
!***********************

!FULL WIDTH HALF MAXIMUM CALCULATION
z=0
y=0
boundInt = x0
WRITE(6,*) 'boundInt:   ', x0

ALLOCATE(gaussian(-boundInt:boundInt)) 
ALLOCATE(fwhm(0:steps)) 

DO k=1, steps
!ASSIGN ONE AXIS TO SEPARATE ARRAY TO DETERMINE MAX VAL
	DO x=-boundInt,boundInt
		gaussian(x)=array(x,y,z,k)
	END DO 

!Determine half maximum width
	tempPeak = MAXVAL(gaussian(-boundInt:boundInt))
	tolerance = 2.0

	WRITE(6,*) k,'tempPeak:   ', tempPeak
	
	!iterating from the lowest x-axis value
	i=-boundInt
	DO
		IF (gaussian(i)>tempPeak/2.0 .AND. gaussian(i)<tempPeak/2.0+tolerance) THEN
			lowVal = i
			WRITE(6,*) k,'lowVal found:  ', lowVal
			EXIT
		ENDIF
		IF (i==boundInt) THEN
		!	WRITE(6,*) k,'no FWHM found'
			EXIT
		ENDIF
		i=i+1
	END DO

	!iterating from the highest x-axis value
	i=boundInt
	DO
		IF (gaussian(i)>tempPeak/2.0 .AND. gaussian(i)<tempPeak/2.0+tolerance) THEN
			highVal = i
			WRITE(6,*) k,'highVal found:  ', highVal
			EXIT
		ENDIF
		IF (i==-boundInt) THEN
		!	WRITE(6,*) k,'no FWHM found'
			EXIT
		ENDIF
		i=i-1
	END DO
	
	fwhm(k) = real(highVal-lowVal)*stepDist
END DO

!***********************
WRITE(6,*) 'CHECK   5'
!***********************

!OUTPUT: print all FWHM values to file
DO i = 1, steps
    WRITE(16,*) real(i)*stepTime, fwhm(i) 
END DO
!OUTPUT: all initial & end values (loop cycles over x-position & then y-position)
	DO k=0,1
		DO z=1-z0,z0-1
			DO y=1-y0,y0-1	
				DO x=1-x0,x0-1	
				WRITE(12+k,*) xAxis(x),yAxis(y),zAxis(z),array(x,y,z,(k*steps))
				END DO 
				WRITE(12+k,*)
			END DO
			WRITE(12+k,*)
		END DO
		WRITE(12+k,*)
	END DO

!OUTPUT:initial & end values for horizonal slices 3D space (loop cycles over x-position & then y-position)
	stepSize = 4
	DO k=0,1
		DO z=1-z0,z0-1-stepSize, stepSize
			DO y=1-y0,y0-1	
				DO x=1-x0,x0-1
					WRITE(14+k,*) xAxis(x), yAxis(y),zAxis(z),array(x,y,z,(k*steps))
				END DO 
				WRITE(14+k,*)
			END DO
			WRITE(14+k,*)
		END DO
	END DO

!OUTPUT:temperature distribution for a slice of the 3D space, over time (loop creates frames in time - 'samples' the data)
	flag = 3
	DO k=0,flag
		z=0	
		DO y=1-y0,y0-1	
			DO x=1-x0,x0-1
				WRITE(30+k,*) xAxis(x), yAxis(y),zAxis(z),array(x,y,z,int(k*steps/flag))
			END DO 
			WRITE(30+k,*)
		END DO
		WRITE(30+k,*)

		bound = 4.0
		DO x=int(-x0+x0/bound),int(x0-1),int(2*x0/bound)	
			DO z=0,z0-1
				DO y=1-y0,y0-1
					WRITE(30+k,*) xAxis(x), yAxis(y),zAxis(z),array(x,y,z,int(k*steps/flag))
				END DO
				WRITE(30+k,*)
			END DO
			WRITE(30+k,*)
		END DO
		WRITE(30+k,*)
	END DO

!***********************
WRITE(6,*) 'CHECK   6'
!***********************

END PROGRAM tempDiffusion