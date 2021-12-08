all:	
	gfortran -O -c forward_rf_subr.f   
	gfortran -o forward_rf  forward_rf.f  forward_rf_subr.o
	rm *.o  
run:
	forward_rf 
	python writesac.py
	#rm SYNRF_t_L_Q.out
        
        
