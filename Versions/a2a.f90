MODULE randomiser

PRIVATE
PUBLIC :: rand, ran_exp

CONTAINS
	REAL FUNCTION rand(iseed,first)
		!
		!  This function returns a pseudo-random number for each invocation.
		!  It is an f90 adaptation of an
		!  FORTRAN 77 adaptation 
		!  by Dick Valent and Fred Clare
		!  http://www.cisl.ucar.edu/zine/96/spring/articles/3.random-6.html
		!  of the "Integer Version 2" minimal 
		!  standard number generator whose Pascal code appears in the article:
		!
		!     Park, Steven K. and Miller, Keith W., "Random Number Generators: 
		!     Good Ones are Hard to Find", Communications of the ACM, 
		!     October, 1988.
		!
		  implicit none
		  integer, parameter :: MPLIER=16807
		  integer, parameter :: MODLUS=2147483647
		  integer, parameter :: MOBYMP=127773
		  integer, parameter :: MOMDMP=2836
		  integer :: hvlue,lvlue,testv,nextn,first,iseed
		  save nextn
		!
		  if(first == 0) THEN
		    nextn=iseed
		    first=1
		  endif
		!
		  hvlue=nextn/mobymp
		  lvlue=mod(nextn,mobymp)
		  testv=mplier*lvlue-momdmp*hvlue
		  if(testv > 0)then
		    nextn=testv
		  else
		    nextn=testv+modlus
		  endif
		  rand = real(nextn)/real(modlus)
	END FUNCTION rand

	REAL FUNCTION ran_exp(iseed,first,lambda)
		

		  implicit none
		  integer :: iseed,first
		  real :: dum,lambda, rand(iseed,first)

		1 continue
		  dum=rand(iseed,first)
		  if(dum == 0.0d00)goto 1
		  ran_exp=-log(dum)*lambda
	END FUNCTION ran_exp
END MODULE randomiser

MODULE alpha
CONTAINS
	SUBROUTINE pi_est(iseed)
				
		USE randomiser

		IMPLICIT NONE
		INTEGER :: n_circ, n, i, iseed, first
		REAL :: x, y, area, pi, error
		
			!format of file
			OPEN(UNIT=20, FILE="euler.dat")
			WRITE(20,'(a9, a11, a10)') '# Time    ', '# Position', '# Velocity'	
			WRITE(20,'(a9, a11, a10)') '# t (s)   ', '# x (m)   ', '# v (m/s) '

			!RAND NUMBER SEQ. GENERATOR
	
			!set number of iterations
			n = 10000
			first = 0		
			iseed = 894989
	
			DO i=1, n
				!set random number to desired range
				x = 2.0*(rand(iseed,first)-0.5)		
				y = 2.0*(rand(iseed,first)-0.5)
				!check if value lies inside circle, add to counter n_circ
				IF(x**2.0+y**2.0<1.0) n_circ = n_circ + 1
				!END IF
			END DO
			!calculate fractional area, hence estimate of area
			area = real(n_circ)/real(n)
			pi = 22.0/7.0
			error = area / pi	
	
			WRITE(6,*) 'area    ', "pi      ", "error   "
			WRITE(6,*) area, pi, error
	END SUBROUTINE
	

	!random distance generator
!std dev of distance
	SUBROUTINE d_path(iseed)

		USE randomiser

		IMPLICIT NONE
		INTEGER :: n_scatt, i, j, n, iseed, first
		REAL :: dist, lambda, std_dv, avg_d, avg_d2, tot_d, tot_d2, err_scatt, ratio
			

			OPEN(UNIT=20, FILE="path.dat")

			!declaration of all variables	
			n = 1000
			j = 100
			n_scatt = 100
			dist = 0.0
			lambda = 1000.0
			tot_d = 0.0
			avg_d = 0.0
			tot_d2 = 0.0 
			
			!loop for changing: n_scatt
			DO n_scatt=1,1000
			
				!loop for all particles (n)
				DO i=1,n
				!loop for 1 particle all events
					dist = 0.0
				
					DO j=1, n_scatt
						dist = dist+ran_exp(iseed, first, lambda)
					END DO
					tot_d = tot_d + dist
					tot_d2 = tot_d2 + dist**2.0
				END DO	
				
				!take avg dist per particle
				avg_d = tot_d/real(n)
				avg_d2 = tot_d2/real(n)
				!std dev per particle
				std_dv = (abs(-avg_d**2.0+avg_d2))**0.5
				!finding congruent error values
				err_scatt = (real(j)*lambda)/(real(i*j)**0.5)

				!find ratio to determine correlation 
				ratio = std_dv/avg_d
				!write ratio vs scattering events to file
				WRITE(20,*) log(ratio), log(real(n_scatt))
			END DO 
	
			!declare values obtained
			WRITE(6,*) 'Standard Deviation:', std_dv
			WRITE(6,*) 'Average Particle Path:', avg_d			
			WRITE(6,*) 'MC Error:', err_scatt
			WRITE(6,*) 'File output on: path.dat'

	END SUBROUTINE d_path
END MODULE alpha


PROGRAM task_1


USE alpha

INTEGER :: iseed



iseed = 894989


CALL pi_est(iseed)
CALL d_path(iseed)



	
END PROGRAM









