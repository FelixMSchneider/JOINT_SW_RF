
JOINT_RF_SW
===========

This code is from the Joint-Inversion Code for Receiver function and Surface wave 

Herrmann, R. B. (2013) Computer programs in seismology: 
An evolving tool for instruction and research, 
Seism. Res. Lettr. 84, 1081-1088, doi:10.1785/0220110096 


I extracted the relevant programs (from VOLIV + shwmod96 from VOLV) for this purpose and added some small features.

Features include:

* Support for L-Q-T Receiver functions
* Normalization of Q-Receiver function with height of L-spike
* Add additional High-Pass Butterworth filter
* simple forward routines (L-Q-T Receiver functions and Dispersion curves) are extracted from the code

The High-Pass Butterwoth filter is controled by USER5 header value of the input SAC Files
For L-Q-T support a near surface Vs must be provided in USER6 header value of the input SAC Files

as in the Original Version of the code

USER4 and B values of the SAC-header control the gauss-parameter and the delay (time before P-onset) of the Receiver function

This slim Version does not need X11, but python and obspy should be installed in order to use the scripts for preparing data.



