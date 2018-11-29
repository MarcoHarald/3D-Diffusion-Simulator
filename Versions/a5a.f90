PROGRAM Bae
        IMPLICIT NONE
	INTEGER, PARAMETER :: m = 1050
	REAL :: prior(0 : m), p(0 : m), tau_arr(0 : m), decay_meas_t(1:10)

	INTEGER :: i, J_MEAS
	REAL :: dtau, norm_const, dx

	!SETTING INITIAL VALUES
	dx = 1.0/(100.0-0.0)       ! 1/(tMAX-tMIN)
	dtau= 0.1 	
	
	!INITIALISE VALUE OF PRIOR
	prior = 0.0
	DO i=10, 1000
		prior(i) = 1.0/(100.0-1.0)
	END DO
        
	!CONNECT FILES FOR INPUT/OUTPUT
	OPEN(UNIT=11, FILE="bayes_data.txt")
	OPEN(UNIT=12, FILE="bayes_data_out.txt")
        
        !SET ACTUAL TIME VALUES FOR TIME VARIABLE
	DO i=0, m
		tau_arr(i)=real(i)*dtau
	END DO

	!READ VALUES FROM INPUT FILE
	DO i=1, 10
		READ(11,*) decay_meas_t(i)
	END DO
	
	!INITIALISE VALUE OF P (PROBABILITY)
	p = prior	

	!IMPROVE VALUES USING BAYES THEOREM  
	DO j_meas=1, 10
		DO i=1, m	
			p(i)=p(i)*exp(-decay_meas_t(j_meas)/tau_arr(i))/tau_arr(i)	
		END DO
	END DO 

	!NORMALISATION OF VALUES
	norm_const = (sum(p)-0.5*(p(0)+p(m)))*dtau
	p = p/norm_const

	!WRITE ALL VALUES TO FILE
        do i=0,m
        write(12,*) tau_arr(i), prior(i), p(i)
        end do

	!CLOSE FILE OUTPUT TO PREVENT ERROR
        close(11)
        close(12)


END PROGRAM
