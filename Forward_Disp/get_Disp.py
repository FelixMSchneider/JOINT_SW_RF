import os

os.system("call_surf > surf.out")


f=open("surf.out")
A=f.readlines()
f.close()

os.system("rm surf.out")


Tg=[float(a.split()[1]) for a in A if a.split()[0]=="GROUP"]
Vg=[float(a.split()[2]) for a in A if a.split()[0]=="GROUP"]

Tp=[float(a.split()[1]) for a in A if a.split()[0]=="PHASE"]
Vp=[float(a.split()[2]) for a in A if a.split()[0]=="PHASE"]



f=open("surf_disp", "w")


for i,T in enumerate(Tp):
    print("SURF96 R C T   0      %7.4f     %7.4f     0.0000" % (T,Vp[i]), file=f)
for i,T in enumerate(Tg):
    print("SURF96 R U T   0      %7.4f     %7.4f     0.0000" % (T,Vg[i]), file=f)



f.close()
