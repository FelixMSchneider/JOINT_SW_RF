import numpy as np
import pylab as plt
import obspy
from obspy.core.util.attribdict import AttribDict
import os

def Gauss(freq,a):
     import numpy as np
     omega=np.pi*2*freq
     return np.exp(-0.25*(omega/a)**2)

def gausfil(data,delta,a):
    from scipy.fftpack import irfft,rfft,rfftfreq
    freq=rfftfreq(len(data), delta)
    gi=Gauss(freq,a)
    fftdata=rfft(data)
    return irfft(gi.T*fftdata)


f=open("INPUT_PARAMETER.dat")
A=f.readlines()
thp=[float(a.split()[0]) for a in A if a.split()[2]=="thp"][0]
rayp=[float(a.split()[0]) for a in A if a.split()[2]=="rayp"][0]
svs=[float(a.split()[0]) for a in A if a.split()[2]=="svs"][0]



plot=True
outputfile=None
outfolder="./"
outfolder="./SAC"
outputfile="RFs.png"
gaus_parameters=[0.25,0.5,1.0,1.5,3.0]

os.system("mkdir -p "+outfolder)

f=open("SYNRF_t_L_Q.out", "r")
data=f.readlines()
f.close()
times   = np.array([float(d.split()[0]) for d in data])
LRFdata = np.array([float(d.split()[1]) for d in data])
QRFdata = np.array([float(d.split()[2]) for d in data])

delta=times[1]-times[0]
delay=10.0


tr=obspy.Trace()
tr.data=QRFdata
    
    
#set headers
tr.stats.delta=delta
tr.stats.station="SYN"
tr.stats.sac=AttribDict()
tr.stats.sac.user4=rayp
tr.stats.sac.user5=thp
tr.stats.sac.user6=svs
tr.stats.sac.b=-delay




for gauspar in gaus_parameters:

    trf=tr.copy()
    GLRFdata=gausfil(LRFdata,delta,gauspar)
    norm=GLRFdata.max()
    GLRFdata/=norm
    trf.data=gausfil(tr.data,delta,gauspar)/norm
    trf.stats.sac.user0=gauspar
    trf.write(outfolder + "/" + "QRF_"+str(gauspar)+".SAC", format="SAC")

    if plot:
        plt.plot(trf.times()-delay, GLRFdata, "b-")
        plt.plot(trf.times()-delay, trf.data, "k-")
    
    
if plot:
    plt.ylim(-0.2,1.1)
#    plt.xlim(-10,40)

if outputfile == None: 
    plt.show()
else:
    plt.savefig(outputfile)

