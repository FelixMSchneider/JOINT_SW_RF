
JOINT_RF_SW
===========

This code is taken and modified from Computer programs in seismology (Rob Hermann) 

&nbsp; Herrmann, R. B. (2013):  
&nbsp; *Computer programs in seismology: An evolving tool for instruction and research,*    
&nbsp; Seism. Res. Lettr. 84, 1081-1088, doi:10.1785/0220110096


I extracted the relevant programs (from VOLIV + shwmod96 from VOLV) for this purpose and added some small features.
This slim Version does not need X11, but python and obspy should be installed in order to use the scripts for preparing data.

Features include:

* Support for L-Q-T Receiver functions
* Normalization of Q-Receiver function with height of L-spike
* Add additional High-Pass Butterworth filter
* simple forward routines (L-Q-T Receiver functions and Dispersion curves) are extracted from the code

The High-Pass Butterwoth filter is controled by USER5 header value of the input SAC Files. <br>
For L-Q-T support a near surface Vs must be provided in USER6 header value of the input SAC Files.

As in the Original Version of the code

USER4 and B values of the SAC-header control the gauss-parameter and the delay (time before P-onset) of the Receiver function




