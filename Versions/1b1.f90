PROGRAM pig
CHARACTER :: word, term

DO i=1,10
	term = char(i)
	word = "endVal" // term // ".dat"
	WRITE(6,*) word
END DO

END PROGRAM

