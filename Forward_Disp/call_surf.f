        program call_surf


c
c   surfdisp96(thkm,vpm,vsm,rhom,nlayer,iflsph,iwave,
c     &                       mode,igr,kmax,t,cg)
c
c

        implicit none
        character mname*80
        character title*80
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs

        integer mmax,i,idimen,icnvel,ierr,iunit,iiso,iflsph, kmax

        common/modlly/mmax


        integer NP
        parameter (NP=60)
        double precision t(NP),cg(NP)



c set parameters

        mname='INPUT_MODEL'   

c       read model


        call getmod(2,mname,mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)

        write(*,*) title 

c-----
c       make sure that we use 1/Q
c-----

       do i=1,mmax
            if(qa(i).lt.0.0)qa(i) = 0.0
            if(qb(i).lt.0.0)qb(i) = 0.0
            if(qa(i) .gt. 1.0)qa(i) = 1.0/qa(i)
            if(qb(i) .gt. 1.0)qb(i) = 1.0/qb(i)
            if(frefp(i) .le. 0.0)frefp(i) = 1.0
            if(frefs(i) .le. 0.0)frefs(i) = 1.0
        enddo


c----- parameters
c     thkm, vpm, vsm, rhom: model for dispersion calculation
c     nlayer - I4: number of layers in the model
c     iflsph - I4: 0 flat earth model, 1 spherical earth model
c     iwave - I4: 1 Love wave, 2 Rayleigh wave
c     mode - I4: ith mode of surface wave, 1 fundamental, 2 first higher, ....
c     igr - I4: 0 phase velocity, > 0 group velocity
c     kmax - I4: number of periods (t) for dispersion calculation
c     t - period vector (t(NP))
c     cg - output phase or group velocities (vector,cg(NP))
c----- 


        kmax=37
        t(1)=5
        t(2)=6
        t(3)=7
        t(4)=8
        t(5)=9
        t(6)=10
        t(7)=11
        t(8)=12
        t(9)=13
        t(10)=14
        t(11)=15
        t(12)=16
        t(13)=17
        t(14)=18
        t(15)=19
        t(16)=20
        t(17)=22
        t(18)=24
        t(19)=26
        t(20)=28
        t(21)=30
        t(22)=32
        t(23)=34
        t(24)=36
        t(25)=38
        t(26)=40
        t(27)=45
        t(28)=50
        t(29)=55
        t(30)=60
        t(31)=65
        t(32)=70
        t(33)=75
        t(34)=80
        t(35)=85
        t(36)=90
        t(37)=95

        call surfdisp96(d,a,b,rho,mmax,iflsph,2,
     &                       1,1,kmax,t,cg)

        do i=1,kmax
           write(*,*) "GROUP", t(i), cg(i)
        enddo
        call surfdisp96(d,a,b,rho,mmax,iflsph,2,
     &                       1,0,kmax,t,cg)

        do i=1,kmax
           write(*,*) "PHASE", t(i), cg(i)
        enddo



        end



