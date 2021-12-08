        program forward_rf
        implicit none
        integer i, isign, k, j, kl, n21, nd2,lw

        integer iunit,iiso, iflsph,idimen,icnvel,ierr
        real rayp, norm, svs
        integer mmx

        real delay, dt, gaussalp, invwgt

        integer n, invdep, iout
        integer nzyear, nzjday, nzhour, nzmin, nzmon, nzday
        logical dop,  outbin

        real freq, df, dfac, fac, fac2
        real wc, hpbutt, thp


c model information
        character mname*80
        character title*80
        character str*80
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs

        integer mmax
        common/modlly/mmax

        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar

        common/depref/refdep
        real refdep



c time series information
        integer NSAMP, NFREQ
        parameter (NSAMP=8192,NFREQ=4097)
        real eata(NSAMP)
        real rfdata(NSAMP, NL+1)
        real qdata(NSAMP, NL+1)
        real ldata(NSAMP, NL+1)


        complex Q(NL+1), Zin
        complex L(NL+1)


        real rdin, cosi, sini, rot11, rot12, rot21, rot22

       common/damp/alpha,ieqex
            real alpha
            integer ieqex

c set parameters

         open(5,file='INPUT_PARAMETER.dat', status='old')

         read(5,'(a)') str
         call lofw(str,lw)
         write(mname,'("./",a,a1)') str(1:lw),char(0)


        read(5,*) rayp
        read(5,*) svs
        read(5,*) delay
        read(5,*) dt
        read(5,*) gaussalp
        read(5,*) n
        read(5,*) thp

        if (thp.gt.0) then
            wc= (1.0/thp) * 6.2831853 
        else
            wc=0.0
        endif


        dop = .true.
        
        outbin = .true.




c       read model


        call getmod(2,mname,mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)

        write(*,*) "" 
        write(*,*) "compute Q-RF and L-RF using input model-file ",mname
        write(*,*) ""

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
c-----
c       define the deepest layer for which a real P-ray is possible
c-----
        call getmmxrayp(mmx,rayp)

c-----
c       ensure that the number of points is a power of 2
c-----

        call npow2(n)

 
c----- 
c       generate a zero phase pulse, with a zero at the  
c       Nyquist Frequency 
c----- 
        delay = abs(delay) 
c----- 
c       define the alpha parameter for complex frequency 
c----- 
        alpha = 2.5/(n*dt)
c-----
c       note this will not work for symmetric pulses with 0 lag ???
c       we remove the lag at a latter stage
c-----

c FS_comment: create initial pulse


c-----
c       now do all the work
c-----
c-----
c       now process in the frequency domain
c-----
        df = 1.0/(n*dt)
        n21 = n / 2 + 1
        do 4000 i=1,n21
            freq=(i-1)*df
            Zin = cmplx(1.0,0.0)

            iout = 4
            invdep=1
            call excit(freq,dop,rayp,Q,Zin,iout,invdep,mmx,svs)
c
c           calculate L-component spike directly, which is a spike at 0
c           (needed for normalization of Q)
c
            do k=1,mmax+1
               L(k) =  dcmplx(1.0d+00, 0.0d+00)
            enddo


            do 4100 kl=1,mmax+1
            fac = - 6.2831853*freq*(delay-dt)
            Q(kl) = Q(kl) * cmplx(cos(fac),sin(fac))
            L(kl) = L(kl) * cmplx(cos(fac),sin(fac))
            if(i.eq.1)then
                qdata(1,kl) = real(Q(kl))
                ldata(1,kl) = real(L(kl))
            else if(i.eq.n21)then
                qdata(2,kl) = real(Q(kl))
                ldata(2,kl) = real(L(kl))
            else
                j = 2*i -1
                k = j + 1
                qdata(j,kl) =  real(Q(kl))
                qdata(k,kl) = aimag(Q(kl))
                ldata(j,kl) =  real(L(kl))
                ldata(k,kl) = aimag(L(kl))
            endif

 4100       continue
 4000   continue

c-----
c       now inverse Fourier Transform everything
c       this is not too bad since we can use Fortran array indexing
c----- 
        do 4200 kl=1,mmax+1
c-----
c           inverse Fourier transform
c-----
            isign = +1
            nd2 = n / 2
            call realft(qdata(1,kl),eata,nd2,isign,dt)
            call realft(ldata(1,kl),eata,nd2,isign,dt)
c-----
c           undamp the time series
c-----
            fac = exp(-alpha*delay)
            dfac = exp(alpha*dt)
            do 425 i = 1,n
                qdata(i,kl)= qdata(i,kl) * fac
                ldata(i,kl)= ldata(i,kl) * fac
                fac = fac * dfac
  425       continue



c-----
c           now Gaussian filter and HP
c-----
            isign = -1
            df = 1.0/(n*dt)
            nd2 = n / 2


c           Fourier transform again
            call realft(qdata(1,kl),eata,nd2,isign,dt)
            call realft(ldata(1,kl),eata,nd2,isign,dt)
            do 426 i=1,n21
                freq=(i-1)*df
                fac = (6.2831853*freq)/(2.0*gaussalp)
                if(fac.gt.25.0)then
                    fac = 0.0
                else
                    fac = exp( - fac * fac)
                endif
                   
                if(wc.gt.0.0) then
                hpbutt=0.0
                else
                hpbutt=1.0
                endif

                if(freq.gt.0) then
                hpbutt = 1/(1+(wc/(6.2831853*freq))**(2*3))
                endif

                fac2=fac*hpbutt

c           apply hp butterworth filter only on q-component
c-----
c           multiplication since we
c           multiply a complex number by a real
c-----
                if(i.eq.1)then
                    qdata(1,kl) = qdata(1,kl) * fac2
                    ldata(1,kl) = ldata(1,kl) * fac
                else if(i.eq.n21)then
c-----
c               Source pulse has Zero at Nyquist
c-----
                    qdata(2,kl) = qdata(2,kl) * fac2
                    ldata(2,kl) = ldata(2,kl) * fac
                else
                    j = 2*i -1
                    k = j + 1
                    qdata(j,kl) = fac2 * qdata(j,kl)
                    qdata(k,kl) = fac2 * qdata(k,kl)
                    ldata(j,kl) = fac * ldata(j,kl)
                    ldata(k,kl) = fac * ldata(k,kl)
                endif
  426       continue
c-----
c           inverse Fourier transform
c-----
            isign = +1
            nd2 = n / 2
            call realft(qdata(1,kl),eata,nd2,isign,dt)
            call realft(ldata(1,kl),eata,nd2,isign,dt)
 4200   continue

c
c         normalize q-component by maximum of l-component
c
           norm=0.0
           do 427 i = 1,n
                if(ldata(i,mmax+1)>norm) then
                   norm=real(ldata(i,mmax+1))
                endif
  427       continue
           


           do 429 kl=1,mmax+1
           do 430 i = 1,n
           ldata(i,kl)=ldata(i,kl)/norm
           qdata(i,kl)=qdata(i,kl)/norm
  430       continue
  429      continue



            open(9,file='SYNRF_t_L_Q.out',status='unknown')


       do 5000 i=1,n
       write (9,*) i*dt-delay, ldata(i,mmax+1), 
     1              qdata(i,mmax+1)
 5000 continue


 


        end



