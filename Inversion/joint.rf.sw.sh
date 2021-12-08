#!/bin/sh

TAIL=`tail -1 rftn.lst`
sta=`saclst kstnm f $TAIL|awk '{print $2}'`

BIN=../bin/

workdir=`pwd`
Dispdir=Disp
Rfdir=Rf

#set vpvs of crust
vpvs=1.8
vp=`echo 4.48*$vpvs|bc -l`

sed 's/X_VP_X/'$vp'/' Line.mod > start.mod


OVeldir=${workdir}/Out/Vel
OResdir=${workdir}/Out/Res
ODispdir=${workdir}/Out/Disp
ORfdir=${workdir}/Out/Rf

mkdir -p ${OVeldir}
mkdir -p ${ODispdir}
mkdir -p ${ORfdir}
mkdir -p ${OResdir}

#####
#	clean up
#####
    $BIN/joint96 39
    $BIN/joint96 36 1
#####
#	set the time window for the RFTN inversion
#	to -5, 20  The 20 is later than the bounce
#	The begin time of -5 is OK since we
#	only use the GAUSS ALPHA=1.0 and 2.5 . If
#	the 0.5 were used we may have to start earlier
#	because of the longer width of the pulse ventered at zero lag
#####
	$BIN/joint96 33 -10
	$BIN/joint96 34 40
#####
#	Std error of fit floor for RFTN
#####
    $BIN/joint96 42 0.0005
#####
#	Std error of fit floor for velocity disp
#####
    $BIN/joint96 40 0.05


# give all filters the same weight
$BIN/joint96 50 1 1
$BIN/joint96 50 2 1
$BIN/joint96 50 3 1
$BIN/joint96 50 4 1
$BIN/joint96 50 5 1




#####
#	Layer weighting
#####
	$BIN/joint96 31 68 0.8
	$BIN/joint96 31 69 0.6
	$BIN/joint96 31 70 0.4
	$BIN/joint96 31 71 0.2
	$BIN/joint96 31 72 0.1

#####
#	also smooth the upper crust a bit
#####
	$BIN/joint96 31  1 0.5
	$BIN/joint96 31  2 0.5
	$BIN/joint96 31  3 0.6
	$BIN/joint96 31  4 0.7
	$BIN/joint96 31  5 0.8
	$BIN/joint96 31  6 0.9
#####
#	set the p factor (joint weighting between the two data sets)
#####
	$BIN/joint96 43 1
#####
#	start the first inversion with a slightly higher damping
#	do avoid an overshoot in the first model estimate
#####
	$BIN/joint96 32 10
	$BIN/joint96 37 2 1 2 6
#####
#	set the p factor
#	do more inversions
#####
	$BIN/joint96 43 0.9
	$BIN/joint96 32 5
	$BIN/joint96 37 2 1 2 6

#####
#	set the p factor
#	do more inversions
#####
	$BIN/joint96 43 0.8
	$BIN/joint96 32 1
	$BIN/joint96 37 2 1 2 6
#####
#	set the p factor
#	do more inversions
#####
	$BIN/joint96 43 0.7
	$BIN/joint96 32 0.5
	$BIN/joint96 37 4 1 2 6
#####
#	set the p factor
#	do 10 more inversions
#####
	$BIN/joint96 43 0.5
	$BIN/joint96 32 0.1
	$BIN/joint96 37 30 1 2 6

#####
#	get the current model
#####
    $BIN/joint96 1 2 28 end.mod
    cp end.mod $OVeldir/Ve.${sta}.end.mod

    $BIN/joint96 1 2 27 dsp.dat
    mv dsp.dat $ODispdir/Dp.${sta}.dat

    $BIN/joint96 1 2 29 res.mod
    mv res.mod $OResdir/Re.${sta}.mod
#####
#	plot up the receiver functions
#####
    $BIN/rftnpv96
    $BIN/plotnps -F7 -W10 -EPS -K < RFTNPV96.PLT > $ORfdir/RFTN.${sta}.ps
	ps2raster $ORfdir/RFTN.${sta}.ps -Tj -A -P
#####
#	plot up the dispersion
#####
    $BIN/srfphv96
    $BIN/plotnps -F7 -W10 -EPS -K < SRFPHV96.PLT > $ODispdir/Disp.${sta}.ps
	ps2raster $ODispdir/Disp.${sta}.ps -Tj -A -P
#####
#	plot up the resolution kernel
#####
    $BIN/srfphr96
    $BIN/plotnps -F7 -W10 -EPS -K < SRFPHR96.PLT > $OResdir/Resolution.${sta}.ps
	ps2raster $-F7 -W10 -EPS -K < SRFPHR96.PLT-ZMAX 200 -K -1 -LEG Line.mod tmpmod96.??? end.modOResdir/Resolution.${sta}.ps -Tj -A -P
#####
#	compare the individual models from the inversion
#	to the end model
#####
    $BIN/shwmod96 -ZMAX 200 -K -1 -LEG Line.mod tmpmod96.??? end.mod
    $BIN/plotnps -F7 -W10 -EPS -K < SHWMOD96.PLT > $OVeldir/Model.${sta}.ps
	ps2raster $OVeldir/Model.${sta}.ps -Tj -A -P


#####
#	output the parameters of inversion
#####
    $BIN/joint96 45 > Weighting.Parameters
    $BIN/joint96 47 > Inversion.Controls
    $BIN/joint96 49 > Rftn.Information

####
#	modeling rf
####

    rm end.mod
    rm *.tmp
    rm tmp*
    rm *PLT
 
