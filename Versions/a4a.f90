PROGRAM Bae
        IMPLICIT NONE
	INTEGER, PARAMETER :: m = 1050
	REAL :: prior(0 : m), p(0 : m), tau_arr(0 : m), decay_meas_t(1:10)

	INTEGER :: i, J_MEAS
	REAL :: dtau, norm_const

	!SETTING INITIAL VALUES
	dx = 1.0/(100.0-0.0)       ! 1/(tMAX-tMIN)
	dtau= 0.1 	
	
	prior = 0.0
	DO i=10, 1000
		prior(i) = 1.0/(100.0-1.0)
	END DO
        
	


	OPEN(UNIT=11, FILE="bayes_data.txt")
	OPEN(UNIT=12, FILE="bayes_data_out.txt")
        
        
	DO i=0, m
		tau_arr(i)=real(i)*dtau
	END DO
        
        do i=0,m
        write(12,*) tau_arr(i), prior(i)
        end do

	!READ VALUES FROM INPUT FILE
	DO i=1, 10
		READ(11,*) decay_meas_t(i)
	END DO

	p = prior	
	WRITE(6,*) decay_meas_t

	!IMPROVE VALUES USING BAYES THEOREM  
	DO j_meas=1, 10
		DO i=1, m	
			p(i)=p(i)*exp(-decay_meas_t(j_meas)/tau_arr(i))/tau_arr(i)	
		END DO
	END DO 

 
        !do i=0,m
        !write(12,*) tau_arr(i), p(i)
       ! end do

	WRITE(6,*) 'STOP 1'

	!NORMALISATION OF VALUES
	
	norm_const = (sum(p)-0.5*(p(0)+p(m)))*dtau
	
	WRITE(6,*) 'STOP 2'
	WRITE(6,*) norm_const		

	p = p/norm_const

	WRITE(6,*) 'STOP 3'
        close(11)
        close(12)


END PROGRAM
