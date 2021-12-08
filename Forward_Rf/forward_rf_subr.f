        subroutine lofw(word,low)
c 
        integer k,lw, low
        character*(*) word
c
        lw=len(word)
        k=0
        do i=1,lw
          if (word(i:i).eq.' ') go to 99
          k=k+1
        end do
99      low=k
c
        return
        end

        subroutine npow2(npts)
c-----
c       Given npts, determine the N=2**m such that N >= npts
c       return the new ntps
c-----
        integer*4 nsamp, npts
        nsamp = npts
        npts = 1
 1000   continue
        if(npts.ge.nsamp)return
        npts = 2*npts
        go to 1000
        end

        subroutine aten(om,qa,qb,xka,xkb,alpha,a,b,atna,atnb,iwat,
     2      frefp,frefs)
c-----
c       make velocities complex, using Futterman causality operator
c-----
        real*4 qa,qb,alpha,a,b
        complex*16 om,at,atna,atnb,xka,xkb
        real*8 pi, om1p, om1s, oml, fac, pi2
        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar
        real*8 CDABS
        complex*16 CDLOG
c-----
c       reference frequency is fref hz
c-----
        om1p=6.2831853*frefp
        om1s=6.2831853*frefs
        pi2 = 1.5707963
        pi=3.1415927d+00
        if(dokjar)then
c-----
c       Kjartansson Constant Q, causal Q operator
c       Kjartansson, E. (1979). 
c           Constant Q-wave propagation and attenuation,
c       J. Geophys. Res. 84, 4737-4748.
c-----
            gama = atan(qa)/pi
            gamb = atan(qb)/pi
            if(gama.le.0.0)then
                atna = cmplx(1.0,0.0)
            else
                fac = pi2*gama
                rfac = sin(fac)/cos(fac)
                atna = dcmplx(1.0d+00,0.0d+00)/
     1              (( (om/om1p)**dble(-gama) ) *
     2               dcmplx(1.0d+00,-dble(rfac)))
            endif
            if(b.gt.1.0e-04*a)then
                if(gamb.le.0.0)then
                    atnb = cmplx(1.0,0.0)
                else
                    fac = pi2*gamb
                    rfac = sin(fac)/cos(fac)
                    atnb = dcmplx(1.0d+00,0.0d+00)/
     1              (( (om/om1s)**dble(-gamb) ) *
     2               dcmplx(1.0d+00,-dble(rfac)))
                endif
            endif
        else
c-----
c       Futterman Causal Q
c-----
c           low frequency cutoff is 0.01 hz
c-----
            oml=0.062831853d+00
            atna=dcmplx(1.0d+00,0.0d+00)
            atnb=dcmplx(1.0d+00,0.0d+00)
            if(qa.gt.0.0)then
              at=dcmplx(0.0d+00,0.0d+00)
              if(CDABS(om).gt.oml) at=CDLOG(om/om1p)/pi
              if(CDABS(om).le.oml) then
                fac=dsqrt(oml*oml + dble(alpha*alpha))/oml
              at=CDLOG(dcmplx(dble(oml),-dble(alpha))/(om1p*fac))/pi
              endif
              atna=(1.+dble(qa)*at+dcmplx(0.0d+00,dble(qa/2.)))
            endif
            if(qb.gt.0.0)then
              at=dcmplx(0.0d+00,0.0d+00)
              if(CDABS(om).gt.oml) at=CDLOG(om/om1s)/pi
              if(CDABS(om).le.oml) then
                fac=dsqrt(oml*oml + dble(alpha*alpha))/oml
              at=CDLOG(dcmplx(dble(oml),-dble(alpha))/(om1s*fac))/pi
              endif
               atnb=(1.+dble(qb)*at+dcmplx(0.0d+00,dble(qb/2.)))
            endif
        endif
        xka=om/(dble(a)*atna)
        if(b.le.1.0e-04*a)then
            iwat = 1
            xkb = dcmplx(0.0d+00,0.0d+00)
        else
            iwat = 0
            xkb=om/(dble(b)*atnb)
        endif
        return
        end

        subroutine hska(aa,w,x,y,z,cosp,cosq,wvno,wvno2,gam,gamm1,rho,
     1      iwat,ex,om2)
        implicit none
c-----
c       command line arguments
c-----
        complex*16 aa(4,4),w,x,y,z,cosp,cosq
        complex*16 wvno,wvno2,gam,gamm1
        real*4 rho
        real*8 ex 
        complex*16 om2
        integer iwat
c-----
c       internal variables
c-----
        complex*16 cpq, gcpq, zw2, gzw2, g1w, g1y, gx
        real*8 dfac
        complex*16 zrho
        integer i,j

        zrho = dcmplx(dble(rho),0.0d+00)
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            do 100 j=1,4
                do 101 i=1,4
                    aa(i,j) = dcmplx(0.0d+00,0.0d+00)
  101           continue
  100       continue
            if(ex.gt.35.0d+00)then
                dfac = 0.0d+00
            else
                dfac = dexp(-ex)
            endif
            aa(1,1) = dfac
            aa(4,4) = dfac
            aa(2,2) = cosp
            aa(3,3) = cosp
            aa(2,3) = -x/zrho
            aa(3,2) = - zrho*w
        else
c-----
c       elastic layer
c-----
c       W = Sa/ra
c       X = ra*Sa
c       Y = Sb/rb
c       Z = rb*Sb
c-----
            cpq = cosp-cosq
            gcpq = gam*cpq
            zw2 = z/wvno2
            gzw2 = gam*zw2
            g1w = gamm1*w
            g1y = gamm1*y
            gx = gam*x
            aa(1,1)=   gcpq + cosq
                aa(1,3)= - wvno * cpq/(zrho*om2)
            aa(1,2)=   wvno*(-g1w+gzw2)
            aa(1,4)=   (wvno2*w-z)/(zrho*om2)
            aa(2,1)=   (gx - wvno2*g1y)/wvno
            aa(2,2)= - gcpq + cosp
            aa(2,3)=   (-x+wvno2*y)/(zrho*om2)
            aa(2,4)= - aa(1,3)
            aa(3,1)=   zrho*om2*gamm1*gcpq/wvno
            aa(3,2)=   zrho*om2*((-gamm1*g1w)+(gam*gzw2))
            aa(3,3)=   aa(2,2)
            aa(3,4)= - aa(1,2)
            aa(4,1)=   zrho*om2*(((gam*gx)/wvno2) - (gamm1*g1y))
            aa(4,2)= - aa(3,1)
            aa(4,3)= - aa(2,1)
            aa(4,4)=   aa(1,1)
        endif
        return
        end

        subroutine var(p,q,ra,rb,w,x,y,z,cosp,cosq,ex,
     1      exa,exl,yl,zl,cosql,iwat)
c     not modified for negative p,q
c     this assumes that real p and real q have same signs
        common/ovrflw/a0,cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz
        complex*16 cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz
        complex*16 p,q,ra,rb,w,x,y,z,cosp,cosq
        complex*16 yl,zl,cosql
        complex*16 eqp,eqm,epp,epm,sinp,sinq
        real *8 a0,pr,pi,qr,qi,fac,qmp,ex,exa,exl
c-----
c       form terms such as cos(p), sin(p), cos(q), sin(q)
c       and cos(p)*cos(q)
c
c       Introduce a factorization of exponentials to
c       make a pseudo floating point ssytem
c
c       ex is the exponent in cosp
c       exl is the exponent in cosq for SH
c       exa is the exponent in cosp*cosq
c-----
        real*8 DREAL
      ex=0.0d+00
      exl = 0.0d+00
      a0=0.0d+00
      pr=dreal(p)
      pi=dimag(p)
      epp=dcmplx(dcos(pi),dsin(pi))/2.
      epm=dconjg(epp)
      ex=pr
      fac=0.0
      if(pr.lt.15.) fac=dexp(-2.*pr)
      cosp=epp + fac*epm
      sinp=epp - fac*epm
      w=sinp/ra
      x=ra*sinp
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            a0 = 1.0d+00
            exa = ex
            cosq = 1.0d+00
            y = 0.0d+00
            z = 0.0d+00
            cosql = 1.0d+00
            yl = 0.0d+00
            zl = 0.0d+00
            exl = 0.0d+00
        else
c-----
c       elastic layer
c-----
            qr=dreal(q)
            qi=dimag(q)
            eqp=dcmplx(dcos(qi),dsin(qi))/2.
            eqm=dconjg(eqp)
            exl=qr
            fac=0.0d+00
            if(qr.lt.15.) fac=dexp(-2.*qr)
            cosql=eqp + fac*eqm
            sinq=eqp - fac*eqm
            yl=sinq/rb
            zl=rb*sinq
c-----
c       form factors for compound P-SV matrix
c-----
            exa=pr + qr
            cpcq=cosp*cosql
            cpy=cosp*yl
            cpz=cosp*zl
            cqw=cosql*w
            cqx=cosql*x
            xy=x*yl
            xz=x*zl
            wy=w*yl
            wz=w*zl
            fac=0.0d+00
            qmp=qr-pr
            if(qmp.gt.-40.) fac=dexp(qmp)
            cosq=cosql*fac
            y=fac*yl
            z=fac*zl
            fac=0.0d+00
            if(exa.lt.60.) a0=dexp(-exa)
        endif
        return
        end

        subroutine excit(freq,dop,rayp,Z,Zin,iout,iptrb,mmx,svs)
c-----
c       Compute the medium response at the given frequency
c
c       freq    R   - frequency in Hz
c       dop     L   - .true.    P - wave incidenc
c                     .false.   S - wave incident
c       rayp    P   - ray parameter in sec/km
c       Z       C   - array of (mmax+1) RFTN's 
c                     mmax + 1 is the unperturbed RFTN
c       Zin     C   - source pulse spectrum at this frequency
c       iout    I   1 output Ur/Uz
c                   2 output Uz
c                   3 output Ur
c       iptrb   I   0 get partials for layer thickness
c                   1 get partials for S velocity
c                   2 get partials for P velocity
c       mmx     I   - deepest layer for which a RFTN is possible, e.g.,
c                     rayp < 1./Vel(P) or c > Vp
c
c-----
        implicit none
c-----
c       command line arguments - define
c-----
        real freq, rayp,svs
        logical dop
        integer iout, iptrb, mmx
        complex Z(*)
        complex Zin
        complex fac
c-----
c       model definition
c-----
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
            real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
            integer mmax
c-----
c       matrix components in layers and boundaries saved
c-----
        common/damp/alpha,ieqex
            real alpha
            integer ieqex

c-----
c       internal variables
c       to do computations we need 1st two columns of 4x4 matrix
c-----
C       complex*16 hv(4,4,NL)
        complex*16 hv(4,2,NL)
        integer k,n
        complex*16 aa(4,4)
        complex*16 paa(4,4)
        real*8 ex(NL)
        real*8 pex(NL)
        real*8 exa, pexa
        complex*16 om, om2, wvno, wvno2
        real omega
        real*8 exl(NL)

        real rdin, cosi, sini, rot11, rot12, rot21, rot22
c-----
c       initialize the solution
c-----
        do k=1,mmax+1
           Z(k) = cmplx(0.0,0.0)
        enddo
c-----
c       if mmx == 0 then there is no ray for this observation
c       and just return so that RFTN and partials are zero
c-----
        if(mmx .eq.0)return
c-----
c       initialize vectors
c       In what follows, k=1,...,mmax-1 refers to the layer k
c                        k=mmax         refers to the halfspace
c                        k=mmax+1       refers to the unperturbed 
c           solution
c       we compute the response for the original and perturbed medium
c       and then compute the partials later by first order one sided
c       differences. Do not worry about correctness of 
c           partial derivative
c       estimate since the estimate will get better as the iterative
c       RFTN inversion converges to a solution
c-----
        do k=1,mmax+1
            hv(1,1,k) = Zin
            hv(1,2,k) = cmplx(0.0,0.0)
            hv(2,1,k) = cmplx(0.0,0.0)
            hv(2,2,k) = Zin
            hv(3,1,k) = cmplx(0.0,0.0)
            hv(3,2,k) = cmplx(0.0,0.0)
            hv(4,1,k) = cmplx(0.0,0.0)
            hv(4,2,k) = cmplx(0.0,0.0)

            exl(k) = 0.0d+00
        enddo
c-----
c       define complex angular frequency and wavenumber
c-----
        omega = 6.2831853*freq
        om =  dcmplx(dble(omega), dble(-alpha))
        om2 = om * om
        wvno  = dcmplx(dble(rayp),0.0d+00)*om
        wvno2 = wvno * wvno



c-----
c       multiply the Haskell matrices from top down
c       taking special care for the perturbed layer
c       n = layer number
c       k = is partial for the layer
c-----
        do n=1,mmx-1
            call gthska (n, ex(n), aa,      om,om2,wvno,wvno2)
            call pgthska(n,pex(n),paa,iptrb,om,om2,wvno,wvno2)

            do k=1,mmax+1
                if(k.eq.n)then
c-----
c                   propagate using perturbed model in layer
c-----
                    call hvm(k,hv,paa)
                    exl(k) = exl(k) + pex(k)
                else
c-----
c                   propagate using original  model in layer
c-----
                    call hvm(k,hv, aa)
                    exl(k) = exl(k) + ex(k)
                endif
            enddo
        enddo
c-----
c       multiply the E sub N sup -1
c-----
            call gteni (mmx, aa,om,om2,wvno,wvno2)
            call pgteni(mmx,paa,om,om2,wvno,wvno2,iptrb)
            do 4000 k=1,mmax+1
                if(k.eq.mmx)then
c-----
c                   propagate using perturbed model in halfspace
c-----
                    call hvm(k,hv,paa)
                else
c-----
c                   propagate using original model in layer
c-----
                    call hvm(k,hv, aa)
                endif
 4000       continue
c-----
c       Now form the solutions - note that I do everything in loops
c       eventually to vectorize
c-----
        if(dop)then
            if(iout.eq.6)then
c-----
c           L-component (Added by Felix Schneider)
c           ROTATION for L and Q component
c           rotation after Svenningsen GJI 2007
c           use s-velocity svs for rotation 
c-----

            rdin=2*asin(rayp*svs)

            cosi=cos(rdin)
            sini=sin(rdin)

            rot11=cosi 
            rot12=sini
            rot21=-sini
            rot22=cosi


c           L = rot11 * Z + rot12 * R
c
c           Z =                       -1 * hv(2,1,k) * fac
c           R = dcmplx(0.0d+00, 1.0d+00) * hv(2,2,k) * fac
c
c           fac = exp(-exl(k))/ (hv(1,1,k)*hv(2,2,k) - hv(2,1,k)*hv(1,2,k)) 

c           L:
            do k=1,mmax+1
              fac=exp(-exl(k))/(hv(1,1,k)*hv(2,2,k)-hv(2,1,k)*hv(1,2,k))
              Z(k) =  ((-1 * rot11 * hv(2,1,k)  +  
     1                  dcmplx(0.0d+00, 1.0d+00) * rot12 * hv(2,2,k)))
     2              * fac
            enddo
            else if(iout.eq.5)then
c-----
c           L-RF L/L (Added by Felix Schneider)
c-----
c           L/L= 1
c           (can be computed directly (without calling of excit function)
            do k=1,mmax+1
               Z(k) =  dcmplx(1.0d+00, 0.0d+00)
            enddo
            else if(iout.eq.4)then
c-----
c           Q-RF Q/L (Added by Felix Schneider)
c           ROTATION for L and Q component
c           rotation after Svenningsen GJI 2007
c           use s-velocity svs for rotation 
c-----

            rdin=2*asin(rayp*svs)

            cosi=cos(rdin)
            sini=sin(rdin)

            rot11=cosi 
            rot12=sini
            rot21=-sini
            rot22=cosi

c           Q= rot21 *Z + rot22 * R 
c           L= rot11 *Z + rot12 * R

c           Z =                       -1 * hv(2,1,k) * fac
c           R = dcmplx(0.0d+00, 1.0d+00) * hv(2,2,k) * fac
c
c           fac = exp(-exl(k))/ (hv(1,1,k)*hv(2,2,k) - hv(2,1,k)*hv(1,2,k)) 
c           (fac is disappears when deviding Q/L)

c           Q/L:
c           rot21 * ( -1* hv(2,1,k)) + rot22 * (dcmplx(0.0d+00, 1.0d+00) * hv(2,2,k)) / 
c           rot11 * ( -1* hv(2,1,k)) + rot12 * (dcmplx(0.0d+00, 1.0d+00) * hv(2,2,k))
c

            do k=1,mmax+1
               Z(k) =  (-1 * rot21 * hv(2,1,k)  +  
     1                  dcmplx(0.0d+00, 1.0d+00) * rot22 * hv(2,2,k))
     2              /  (-1 * rot11 * hv(2,1,k)  +  
     3                   dcmplx(0.0d+00, 1.0d+00) * rot12 * hv(2,2,k))
            enddo
            else if(iout.eq.1)then
c-----
c           UR/UZ
c-----
            do k=1,mmax+1
               Z(k) = - dcmplx(0.0d+00, 1.0d+00) * hv(2,2,k) 
     1             / hv(2,1,k)
            enddo
            else if(iout.eq.2)then
c-----
c           Uz
c-----
            do  k=1,mmax+1
               Z(k) = -  hv(2,1,k)*exp(-exl(k))/
     1             (hv(1,1,k)*hv(2,2,k) - hv(2,1,k)*hv(1,2,k))
            enddo
            else if(iout.eq.3)then
c-----
c           Ur
c-----
            do k=1,mmax+1
               Z(k) = dcmplx(0.0d+00, 1.0d+00) * hv(2,2,k)*exp(-exl(k))/
     1             (hv(1,1,k)*hv(2,2,k) - hv(2,1,k)*hv(1,2,k))
            enddo
            endif
        else
            if(iout.eq.1)then
c-----
c           UR/UZ
c-----
            do k=1,mmax+1
               Z(k) = - dcmplx(0.0d+00, 1.0d+00) * hv(1,2,k) 
     1             / hv(1,1,k)
            enddo
            else if(iout.eq.2)then
c-----
c           UZ
c-----
            do  k=1,mmax+1
            Z(k) = dcmplx(0.0d+00, 1.0d+00) * hv(1,1,k)*exp(-exl(k))/
     1          (hv(1,1,k)*hv(2,2,k) - hv(2,1,k)*hv(1,2,k))
            enddo
            else if(iout.eq.3)then
c-----
c           Ur
c-----
            do  k=1,mmax+1
            Z(k) = -                          hv(1,2,k)*exp(-exl(k))/
     1          (hv(1,1,k)*hv(2,2,k) - hv(2,1,k)*hv(1,2,k))
            enddo
            endif
        endif
c-----
c       do a final cleanup with respect to partials for mmx < mmax
c-----
        do k=1,mmax
           if(k.gt.mmx)then
              Z(k) = cmplx(0.0,0.0)
           endif
        enddo
            
        return
        end

        subroutine gthska(m,ex,aa,om,om2,wvno,wvno2)
        implicit none
        integer m
        real*8 ex
        complex*16 aa(4,4)
        integer iptrb
        complex*16 om, om2, wvno, wvno2

        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
            real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
            integer mmax
        common/damp/alpha,ieqex
            real alpha
            integer ieqex

        complex*16 w,x,y,z,cosp,cosq,gam,gamm1
        complex*16 xka,xkb,ra,rb
        complex*16 atna, atnb,p,q
        complex*16 yl,zl,cosql
        real*8 exa, exb
        integer iwat

        call aten(om,qa(m),qb(m),xka,xkb,
     1      alpha,a(m),b(m),atna,atnb,iwat,
     2      frefp(m),frefs(m))
        gam=dble(b(m))*(wvno/om)
        gam = gam * atnb
        gam = dcmplx(2.0d+00,0.0d+00)*gam*gam
        gamm1 = gam - dcmplx(1.0d+00,0.0d+00)
        ra=CDSQRT(wvno2-xka*xka)
        rb=CDSQRT(wvno2-xkb*xkb)
        p=ra*dble(d(m))
        q=rb*dble(d(m))
        call var(p,q,ra,rb,w,x,y,z,cosp,cosq,
     1          ex,exa,exb,yl,zl,cosql,iwat)
        call hska(aa,w,x,y,z,cosp,cosq,wvno,wvno2,gam,gamm1,rho(m),
     1      iwat,ex,om2)
        return
        end

        subroutine pgthska(n,ex,aa,iptrb,om,om2,wvno,wvno2)
        implicit none
        integer n
        real*8 ex
        complex*16 aa(4,4)
        integer iptrb
        complex*16 om,om2, wvno, wvno2

        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
            real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
            integer mmax
        real hsav, asav, bsav
c-----
c       perturb the model
c-----
        if(iptrb.eq.0)then
            hsav = d(n)
            d(n) = 1.01 * d(n)
        else if(iptrb.eq.1)then
            bsav = b(n)
            b(n) = 1.01 * b(n)
        else if(iptrb.eq.2)then
            asav = a(n)
            a(n) = 1.01 * a(n)
        endif
c-----
        call gthska(n,ex,aa,om,om2,wvno,wvno2)
c-----
c       un perturb the model
c-----
        if(iptrb.eq.0)then
            d(n) = hsav
        else if(iptrb.eq.1)then
            b(n) = bsav
        else if(iptrb.eq.2)then
            a(n) = asav
        endif
        return
        end

        subroutine gteni(m,g,om,om2,wvno,wvno2)
        implicit none
        integer m
        complex*16 g(4,4)
        complex*16 om, om2, wvno, wvno2
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
            real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
            integer mmax
        common/damp/alpha,ieqex
            real alpha
            integer ieqex
c-----
c       get elements of E sub N sup -1 matrix
c-----
        complex*16 xka,xkb,ra,rb,gam,gamm1
        complex*16 atna, atnb
        complex*16 CDSQRT
        integer iwat
c-----
c       set up halfspace conditions
c-----
            call aten(om,qa(m),qb(m),xka,xkb,
     1          alpha,a(m),b(m),atna,atnb,iwat,
     2          frefp(m),frefs(m))
            gam=dble(b(m))*(wvno/om)
            gam = gam * atnb
            gam = dcmplx(2.0d+00,0.0d+00)*gam*gam
            gamm1 = gam - dcmplx(1.0d+00,0.0)
            ra=CDSQRT(wvno2-xka*xka)
            rb=CDSQRT(wvno2-xkb*xkb)
            g(1,1) =   gam/wvno
            g(2,1) = - gamm1/rb
            g(3,1) =   g(1,1)
            g(4,1) =   g(2,1)
            g(1,2) =  - gamm1 / ra
            g(2,2) =   g(1,1)
            g(3,2) =    gamm1 / ra
            g(4,2) = - g(1,2)
            g(1,3) = - 1.0d+00/(rho(m)*om2)
            g(2,3) =   wvno / ( rb * rho(m)*om2)
            g(3,3) =   g(1,3)
            g(4,3) =   g(3,3)
            g(1,4) =   wvno / ( ra * rho(m)*om2)
            g(2,4) =   g(1,3)
            g(3,4) = - wvno / ( ra * rho(m)*om2)
            g(4,4) = - g(2,4)
        return
        end

        subroutine pgteni(n,aa,om,om2,wvno,wvno2,iptrb)
        implicit none
        integer n
        complex*16 aa(4,4)
        complex*16 om, om2, wvno, wvno2
        integer iptrb
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
            real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
            integer mmax
        real hsav, asav, bsav
c-----
c       matrix components in layers and boundaries saved
c-----
c       perturb the model
c       note that if iptrb = 0 for layer thickness we do nothing
c-----
        if(iptrb.eq.1)then
            bsav = b(n)
            b(n) = 1.01 * b(n)
        else if(iptrb.eq.2)then
            asav = a(n)
            a(n) = 1.01 * a(n)
        endif
        call gteni(n,aa,om,om2,wvno,wvno2)
c-----
c       un perturb the model
c-----
        if(iptrb.eq.1)then
            b(n) = bsav
        else if(iptrb.eq.2)then
            a(n) = asav
        endif
        return
        end

        subroutine hvm(n,hv,aa)
c-----
c       FORM hv = aa hv
c-----
        implicit none
        integer NL
        parameter(NL=200)
        complex*16 hv(4,2,NL)
C       complex*16 hv(4,4,NL)
        integer n
        complex*16 aa(4,4)

        complex*16 a11, a21, a31, a41
        complex*16 a12, a22, a32, a42
C       COMPLEX*16 A13, A23, A33, A43
C       COMPLEX*16 A14, A24, A34, A44

        a11 = aa(1,1)*hv(1,1,n) + aa(1,2)*hv(2,1,n) 
     1      + aa(1,3)*hv(3,1,n) + aa(1,4)*hv(4,1,n)
        a12 = aa(1,1)*hv(1,2,n) + aa(1,2)*hv(2,2,n) 
     1      + aa(1,3)*hv(3,2,n) + aa(1,4)*hv(4,2,n)
        a21 = aa(2,1)*hv(1,1,n) + aa(2,2)*hv(2,1,n) 
     1      + aa(2,3)*hv(3,1,n) + aa(2,4)*hv(4,1,n)
        a22 = aa(2,1)*hv(1,2,n) + aa(2,2)*hv(2,2,n) 
     1      + aa(2,3)*hv(3,2,n) + aa(2,4)*hv(4,2,n)
        a31 = aa(3,1)*hv(1,1,n) + aa(3,2)*hv(2,1,n) 
     1      + aa(3,3)*hv(3,1,n) + aa(3,4)*hv(4,1,n)
        a32 = aa(3,1)*hv(1,2,n) + aa(3,2)*hv(2,2,n) 
     1      + aa(3,3)*hv(3,2,n) + aa(3,4)*hv(4,2,n)
        a41 = aa(4,1)*hv(1,1,n) + aa(4,2)*hv(2,1,n) 
     1      + aa(4,3)*hv(3,1,n) + aa(4,4)*hv(4,1,n)
        a42 = aa(4,1)*hv(1,2,n) + aa(4,2)*hv(2,2,n) 
     1      + aa(4,3)*hv(3,2,n) + aa(4,4)*hv(4,2,n)

        hv(1,1,n) = a11
        hv(1,2,n) = a12
        hv(2,1,n) = a21
        hv(2,2,n) = a22
        hv(3,1,n) = a31
        hv(3,2,n) = a32
        hv(4,1,n) = a41
        hv(4,2,n) = a42

        return
        end

        subroutine getmmxrayp(mmx,rayp)
c-----
c       get the maximum layer for which rayp < 1./Vp
c
c       mmx    I   - deepest layer for which rayp < 1./Vp
c                    mmx will be returned in the range [0,mmax]
c                    0 means that no real ray is possible, and thus
c                    the predicted function and the partials will be
c                    set equal to zero.
c                    For mmx >= 1, the partials will be zero for
c                    layers mmx+1, mmax
c       rayp   R   - ray parameter
c-----
        implicit none
        integer mmx
        real rayp
c
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
            real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
            integer mmax

        integer i
c-----
c       safety
c-----
        mmx = 0
c-----
c       do search
c-----
        do i=1,mmax
           if(rayp .lt. 1./a(i))then
               mmx = i
           else
               return
           endif
        enddo
        return
        end
       

        subroutine ftwofft(DATA1, DATA2, FFT1, FFT2, dt, df, N)
c-----
c       E. Oran Brigham
c       The Fast Fourier Transform and Its Applications
c       Prentice Hall, Englewood Cliffs, New Jersey
c       1988
c       ISBN 0-13-307505-2
c-----
c       Find the Discrete Fourier transform of two real time
c       series of length N with just one call to four()
c       Brigham (1988), Sec. 9.3
c-----
c       DATA1   R   - real array 1
c       DATA2   R   - real array 2
c       FFT1    C   - complex transform of DATA1 - output -
c       FFT2    C   - complex transform of DATA2 - output -
c       dt  R   - time sampling interval
c       df  R   - frequency damping interval = 1/(N*dt)
c       N   I   - number of points, power of 2
c-----
        real DATA1(N), DATA2(N)
        complex FFT1(N), FFT2(N)
        real rn, rnmn, in, inmn

        do 1000 i=1,N
            fft1(i) = cmplx( data1(i), data2(i))
 1000   continue
        call zfour(fft1,n,-1,dt,df)
c-----
c       now unwrap
c-----
        n21 = n / 2 + 1
        do 2000 i=1,n21
        if(i.eq.1)then
            rn =  real(fft1(1))
            in = aimag(fft1(1))
            fft1(1) = cmplx(rn, 0.0)
            fft2(1) = cmplx(in, 0.0)
        else
            rn   =  real(fft1(i    ))
            in   = aimag(fft1(i    ))
            rnmn =  real(fft1(n+2-i))
            inmn = aimag(fft1(n+2-i))
            fft1(i) = cmplx(0.5*(rn+rnmn), 0.5*(in-inmn))
            fft2(i) = cmplx(0.5*(in+inmn),-0.5*(rn-rnmn))
            fft1(n+2-i) = conjg(fft1(i))
            fft2(n+2-i) = conjg(fft2(i))
        endif
 2000   continue
        fft1(n21) = cmplx(real(fft1(n21)), 0.0)
        fft2(n21) = cmplx(real(fft2(n21)), 0.0)
        return
        end

        subroutine itwofft(DATA1, DATA2, FFT1, FFT2, dt, df, N)
c-----
c       E. Oran Brigham
c       The Fast Fourier Transform and Its Applications
c       Prentice Hall, Englewood Cliffs, New Jersey
c       1988
c       ISBN 0-13-307505-2
c-----
c       Find the inverse Discrete Fourier transform of two real time
c       series of length N with just one call to four()
c       Brigham (1988), Sec. 9.3
c-----
c       DATA1   R   - real array 1 - output -
c       DATA2   R   - real array 2 - output -
c       FFT1    C   - complex transform of DATA1
c       FFT2    C   - complex transform of DATA2
c       dt  R   - time sampling interval
c       df  R   - frequency damping interval = 1/(N*dt)
c       N   I   - number of points, power of 2
c-----
        real DATA1(N), DATA2(N)
        complex FFT1(N), FFT2(N)
        real rn, rnmn, in, inmn

c-----
c       now wrap
c-----
        n21 = n / 2 + 1
        do 2000 i=1,n
        if(i.eq.1)then
            rn =  real(fft1(1))
            in =  real(fft2(1))
            fft1(1) = cmplx( rn, in)
        else
            in   =  aimag(fft1(i    )) +  real(fft2(i     ))
            inmn = -aimag(fft1(i    )) +  real(fft2(i     ))
            rn   =   real(fft1(i    )) - aimag(fft2(i     ))
            rnmn =   real(fft1(i    )) + aimag(fft2(i     ))
            
            fft1(i)     = cmplx( rn  , in  )
C           fft1(n+2-i) = cmplx( rnmn, inmn)
        endif
 2000   continue

        call zfour(fft1,n,+1,dt,df)
        do 1000 i=1,N
            data1(i) =  real(fft1(i))
            data2(i) = aimag(fft1(i))
 1000   continue
        return
        end

        subroutine realft(data,eata,n,isign,dt)
c-----
c       compute inverse FFT for  2n real FFT with one 
c       N FFT operation isign = +1
c       derived following Brigham (1988), Sec. 9.3
c-----
c       data    R   - array of 2n values
c       n   I   - 2n / 2
c       isign   I   - -1 forward FFT, +1 inverse FFT
c
c       for ISIGN < 0
c       Input :  FORTRAN indexing
c           data(1), data(2), ..., data(2n-1)
c       Output: FORTRAN indexing
c           R(1), R(n21)   where n21 = 2n/2 +1
c           R(2), I(2)
c           ...   ...
c           R(2n),I(2n)
c       we can get others for real time series by Real Even, Imag odd
c
c       for ISIGN > 0
c       Output:  FORTRAN indexing
c           data(1), data(2), ..., data(2n-1)
c       Input : FORTRAN indexing
c           R(1), R(n21)   where n21 = 2n/2 +1
c           R(2), I(2)
c           ...   ...
c           R(2n),I(2n)
c       we can get others for real time series by Real Even, Imag odd
c-----
        real data(n+n)
        real eata(n+n)
        integer n, isign

        real re, ro, ie, io
        real*8 dtheta, s, c, ds, dc

c-----
c       The proper Delta T and Delta F will be imposed later
c       just force the FFT to use ddt = 1 for the FFT, and
c       ddf = 1 for the inverse FFT
c-----
        if(isign.lt.0)then
            ddt = 1.0
            ddf = 0.0
            df = 1.0/(2*n*dt)
            call rfour(data,n,isign,ddt,ddf)
        else
            df = 1.0/(2*n*dt)
            ddt = 0.0
            ddf = 1.0
        endif

c----
c       rearrange
c-----
        dtheta = 3.1415927/real(n)
        dc = cos(dtheta)
        ds = sin(dtheta)
        do 1000 i=1,N

        if(i.eq.1)then
            if(isign.lt.0)then
                eata(1) = (data(1) + data(2))
                eata(2) = (data(1) - data(2))
            else
                eata(1) = (data(1) + data(2)) /2.0
                eata(2) = (data(1) - data(2)) /2.0
            endif
            c = 1.0d+00
            s = 0.0d+00
        else
            jr = 2*i -1
            ji = 2*i
            kr = 2*n + 2*(2 - i) -1
            ki = 2*n + 2*(2 - i) 
            ts = s*dc + c*ds
            tc = c*dc - s*ds
            s = ts
            c = tc

            re =  ( data(jr) +  data(kr) )/2.0
            io =  ( data(ji) -  data(ki) )/2.0
            
            if(isign.lt.0)then
                ro =  ( data(jr) -  data(kr) )/2.0
                ie =  ( data(ji) +  data(ki) )/2.0
                eata(jr) =   re + c*ie - s*ro
                eata(ji) =   io - s*ie - c*ro
                eata(kr) =   re - c*ie + s*ro
                eata(ki) = - io - s*ie - c*ro
            else
                ie = 0.5 * 
     1              ( -s * (data(ji)+data(ki)) + c*(data(jr)-data(kr)))
                ro = 0.5 * 
     1              ( -c * (data(ji)+data(ki)) - s*(data(jr)-data(kr)))
                eata(jr) = re + ro
                eata(ji) = io + ie 
                eata(kr) = re - ro
                eata(ki) = ie - io 
            endif
        endif
 1000   continue
        do i=1,2*n
           if(isign.lt.0)then
               data(i) = eata(i) * dt
           else if(isign.gt.0)then
               data(i) = 2.0 * eata(i) * df
           endif
        enddo
        if(isign.gt.0)then
            call rfour(data,n,isign,ddt,ddf)
        endif
        return
        end

        subroutine rfour(rarr,nn,isign,dt,df) 
c-----
c     THE input is a real array
c     which has numbers stored in memory as
c     R1, I1, R2, I2, ..., Rnn, Inn
c     where nn must be a power of 2 R and I are the real and imaginary
c     parts of the complex number
c
c     For isign -1 this is a complex time series
c     For isign +1 this is a complex frequency series with
c        index 1 (in fortran corresponding to f=0
c              2                              f=df
c            nn/2 + 1                         f = 1/2dt = Nyquist
c            nn - 1                           f = -2df
c             nn                              f = -df
c-----
c     the cooley-tookey fast fourier transform in usasi basic fortran
c     transform(j) = sum(datc(i)*w**((i-1)(j-1)), where i and j run
c     from 1 to nn and w = exp(isign*2*pi*sqrt(-1)/nn).  datc is a one-
c     dimensional complex array (i.e., the real and imaginary parts of
c     datc are located immediately adjacent in storage, such as fortran
c     places them) whose length nn is a power of two.  isign
c     is +1 or -1, giving the sign of the transform.  transform values
c     are returned in array datc, replacing the input datc.  the time is
c     proportional to n*log2(n), rather than the usual n**2
c     rms resolution error being bounded by 6*sqrt(i)*log2(nn)*2**(-b),
c     b is the number of bits in the floating point fraction.
c
c     the program computes df from dt, dt from df and checks to see
c     if they are consistent. In addition, the transforms are multiplied
c     by dt or df to make the results dimensionally correct
c
c     This is a slightly modified version of the original Brenner routine
c     The changes were to add physical dimensions to the transform
c-----
        real rarr(*) 
        n = 2 * nn 
        if(dt.eq.0.0) dt = 1./(nn*df) 
        if(df.eq.0.0) df = 1./(nn*dt) 
        if(dt.ne.(nn*df)) df = 1./(nn*dt) 
        j = 1 
        do  i=1,n,2 
           if(i .lt. j) then
               go to 1
           else
               go to 2
           endif
    1      tempr = rarr(j) 
           tempi = rarr(j+1) 
           rarr(j) = rarr(i) 
           rarr(j+1)=rarr(i+1) 
           rarr(i) = tempr 
           rarr(i+1) = tempi 
    2      m = n/2 
    3      continue
           if(j.le.m) then
               go to 5
           else 
               go to 4
           endif
    4      j = j-m 
           m = m/2 
           if(m.lt.2)then
               go to 5
           else
               go to 3
           endif
    5 j=j+m 
      enddo
        mmax = 2 
C    6 if(mmax-n) 7,10,10 
    6 continue
        if(mmax .lt. n)then
            go to 7
        else if(mmax .ge. n)then
            go to 10
        endif
    7 istep= 2 *mmax 
        theta = 6.283185307/float(isign*mmax) 
        sinth=sin(theta/2.) 
        wstpr=-2.*sinth*sinth 
        wstpi=sin(theta) 
        wr=1.0 
        wi=0.0 
        do  m=1,mmax,2 
           do  i=m,n,istep 
              j=i+mmax 
              tempr=wr*rarr(j)-wi*rarr(j+1) 
              tempi=wr*rarr(j+1)+wi*rarr(j) 
              rarr(j)=rarr(i)-tempr 
              rarr(j+1)=rarr(i+1)-tempi 
              rarr(i)=rarr(i)+tempr 
              rarr(i+1) = rarr(i+1)+tempi 
           enddo
           tempr = wr 
           wr = wr*wstpr-wi*wstpi + wr 
           wi = wi*wstpr+tempr*wstpi + wi 
        enddo
        mmax = istep 
        go to 6 
   10 continue 
        if(isign.gt.0) then
c-----
c           frequency to time domain 
c-----
           do  iiii = 1,n 
               rarr(iiii) = rarr(iiii) * df 
           enddo
        else
c-----
c          time to frequency domain 
c-----
           do  iiii = 1,n 
              rarr(iiii) = rarr(iiii) * dt 
           enddo
        endif
        return 
        end 

        subroutine zfour(zarr,nn,isign,dt,df) 
c-----
c     THE input is a complex array
c     which has numbers stored in memory as
c     R1, I1, R2, I2, ..., Rnn, Inn
c     where nn must be a power of 2 R and I are the real and imaginary
c     parts of the complex number
c
c     For isign -1 this is a complex time series
c     For isign +1 this is a complex frequency series with
c        index 1 (in fortran corresponding to f=0
c              2                              f=df
c            nn/2 + 1                         f = 1/2dt = Nyquist
c            nn - 1                           f = -2df
c             nn                              f = -df

c-----
c     the cooley-tookey fast fourier transform in usasi basic fortran
c     transform(j) = sum(datc(i)*w**((i-1)(j-1)), where i and j run
c     from 1 to nn and w = exp(isign*2*pi*sqrt(-1)/nn).  datc is a one-
c     dimensional complex array (i.e., the real and imaginary parts of
c     datc are located immediately adjacent in storage, such as fortran
c     places them) whose length nn is a power of two.  isign
c     is +1 or -1, giving the sign of the transform.  transform values
c     are returned in array datc, replacing the input datc.  the time is
c     proportional to n*log2(n), rather than the usual n**2
c     rms resolution error being bounded by 6*sqrt(i)*log2(nn)*2**(-b),
c     b is the number of bits in the floating point fraction.
c
c     the program computes df from dt, dt from df and checks to see
c     if they are consistent. In addition, the transforms are multiplied
c     by dt or df to make the results dimensionally correct
c
c     This is a slightly modified version of the original Brenner routine
c     The changes were to add physical dimensions to the transform
c     and to make it all complex
c-----
        complex zarr(*) 
        integer nn, isign
        real dt, df

        complex ztemp
c-----
c       ensure that the dt and df are defined and
c       consistent
c-----
        if(dt.eq.0.0) dt = 1./(nn*df) 
        if(df.eq.0.0) df = 1./(nn*dt) 
        if(dt.ne.(nn*df)) df = 1./(nn*dt) 
c-----
c       now begin the transform
c-----
        jj = 1
        do 5 ii=1,nn 
        if(ii .lt. jj) then
              ztemp = zarr(jj)
              zarr(jj) = zarr(ii)
              zarr(ii) = ztemp
        endif
        mm = nn/2
    3   continue
        if(jj.le.mm) then
            go to 55
        else 
              jj = jj - mm
              mm = mm /2
              if(mm.lt.1)then
                  go to 55
              else
                  go to 3
              endif
        endif
   55   continue
        jj = jj + mm
    5   continue
        mmax = 1 
c-----
    6 continue
        if(mmax .lt. nn)then
            go to 7
        else if(mmax .ge. nn)then
            go to 10
        endif
    7   continue
        istep= 2*mmax 
        theta = 6.283185307/(isign*2.0*mmax) 
        sinth=sin(theta/2.) 
        wstpr=-2.*sinth*sinth 
        wstpi=sin(theta) 
        wr=1.0 
        wi=0.0 
        do 9 mm=1,mmax
              do 8 ii=mm,nn,istep
                    jj=ii+mmax
                    ztemp=cmplx(wr,wi)*zarr(jj)
                    zarr(jj) = zarr(ii) - ztemp
                    zarr(ii) = zarr(ii) + ztemp
    8         continue
c-----
c       use trig relations to compute the next sin/cos
c       without actually calling sin() or cos()
c-----
              tempr = wr 
              wr = wr*wstpr-wi*wstpi + wr 
              wi = wi*wstpr+tempr*wstpi + wi 
    9   continue
        mmax = istep 
        go to 6 
c-----
c       transform is done
c-----
   10   continue 
c-----
c     give the arrays the proper physical dimensions
c-----
        if(isign.lt.0)then
c-----
c             time to frequency domain
c-----
              do  ii = 1,nn
                    zarr(ii) = zarr(ii) * dt
              enddo
        else
c-----
c             frequency to time domain
c-----
              do ii = 1,nn
                    zarr(ii) = zarr(ii) * df
              enddo
        endif
        return
        end
        subroutine getmod(rlun,mname,mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,listmd)
c-----
c       HISTORY
c
c       09 08 2000  gave ierr an initial default value for g77
c       01 13 2001  put in close(lun) if file is not model file
c       03 MAY 2002     Modify to permit read from standard input
c       06 JUL 2005 moved inquire to permit use of STDIN
c
c-----
c       General purpose model input
c       This model specification is designed to be as 
c           general as possible
c
c       Input lines
c       Line 01: MODEL
c       Line 02: Model Name
c       Line 03: ISOTROPIC or ANISOTROPIC or 
c           TRANSVERSELY ANISOTROPIC
c       Line 04: Model Units, First character 
c           is length (k for kilometer
c           second is mass (g for gm/cc), third is time (s for time)
c       Line 05: FLAT EARTH or SPHERICAL EARTH
c       Line 06: 1-D, 2-D or 3-D
c       Line 07: CONSTANT VELOCITY
c       Line 08: open for future use
c       Line 09: open for future use
c       Line 10: open for future use
c       Line 11: open for future use
c       Lines 12-end:   These are specific to the model
c           For ISOTROPIC the entries are
c           Layer Thickness, P-velocity, S-velocity, Density, Qp, Qs,
c           Eta-P, Eta S (Eta is frequency dependence), 
c           FreqRefP, FreqRefP
c-----
cMODEL
cTEST MODEL.01
cISOTROPIC
cKGS
cFLAT EARTH
c1-D
cCONSTANT VELOCITY
cLINE08
cLINE09
cLINE10
cLINE11
c H  VP  VS   RHO   QP  QS   ETAP   ETAS REFP  REFS
c1.0    5.0 3.0 2.5 0.0 0.0 0.0 0.0 1.0 1.0
c2.0    5.1 3.1 2.6 0.0 0.0 0.0 0.0 1.0 1.0
c7.0    6.0 3.5 2.8 0.0 0.0 0.0 0.0 1.0 1.0
c10.0   6.5 3.8 2.9 0.0 0.0 0.0 0.0 1.0 1.0
c20.0   7.0 4.0 3.0 0.0 0.0 0.0 0.0 1.0 1.0
c40.0   8.0 4.7 3.3 0.0 0.0 0.0 0.0 1.0 1.0
c-----
c-----
c       rlun    I*4 - logical unit for reading model file. This
c                 unit is released after the use of this routine
c       mname   C*(*)   - model name - if this is stdin or 
c           STDIN just read
c                 from standard input
c       mmax    I*4 - number of layers in the model, last layer is
c                    halfspace
c       title   C*(*)   - title of the model file
c       iunit   I*4 - 0 Kilometer, Gram, Sec
c       iiso    I*4 - 0 isotropic 
c                 1 transversely anisotropic 
c                 2 general anisotropic 
c       iflsph  I*4 - 0 flat earth model
c                 1 spherical earth model
c       idimen  I*4 - 1 1-D
c               - 2 2-D
c               - 3 3-D
c       icnvel  I*4 - 0 constant velocity
c                 1 variable velocity
c       ierr    I*4 - 0 model file correctly read in
c               - -1 file does not exist
c               - -2 file is not a model file
c                 -3 error in the model file
c       listmd  L   - .true. list the model
c------

        implicit none
        character mname*(*), title*(*)
        integer rlun
        integer*4 mmax, iunit, iiso, iflsph, idimen, icnvel
        integer*4 ierr
        character string*80
        logical listmd
c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)

        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep

        logical ext
        character ftype*80
        integer lun, j, i, irefdp

c-----
c       test to see if the file exists
c-----
        ierr = 0
c-----
c       test for input
c-----
        if(MNAME(1:5).eq.'stdin' .or. mname(1:5).eq.'STDIN')then
c-----
c           do not open anything, use standard output
c-----
            lun = LIN
        else
            lun = rlun
            inquire(file=mname,exist=ext)
            if(.not.ext)then
                ierr = -1
                write(LER,*)'Model file does not exist'
                return
            endif
c-----
c           open the file
c-----
            open(lun,file=mname,status='old',form='formatted',
     1          access='sequential')
            rewind lun
        endif
c-----
c       verify the file type
c-----
c-----
c       LINE 01
c-----
        read(lun,'(a)')ftype
        if(ftype(1:5).ne.'model' .and. ftype(1:5).ne.'MODEL')then
            ierr = -2
            write(LER,*)'Model file is not in model format'
            close(lun)
            return
        endif
c-----
c       LINE 02
c-----
        read(lun,'(a)')title
c-----
c       LINE 03
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'ISO' .or. string(1:3).eq.'iso')then
            iiso = 0
        else if(string(1:3).eq.'TRA' .or. string(1:3).eq.'tra')then
            iiso = 1
        else if(string(1:3).eq.'ANI' .or. string(1:3).eq.'ani')then
            iiso = 2
        endif
c-----
c       LINE 04
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'KGS' .or. string(1:3).eq.'kgs')then
            iunit = 0
        endif
c-----
c       LINE 05
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'FLA' .or. string(1:3).eq.'fla')then
            iflsph = 0
        else if(string(1:3).eq.'SPH' .or. string(1:3).eq.'sph')then
            iflsph = 1
        endif
c-----
c       LINE 06
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'1-d' .or. string(1:3).eq.'1-D')then
            idimen = 1
        else if(string(1:3).eq.'2-d' .or. string(1:3).eq.'2-D')then
            idimen = 2
        else if(string(1:3).eq.'3-d' .or. string(1:3).eq.'3-D')then
            idimen = 3
        endif
c-----
c       LINE 07
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'CON' .or. string(1:3).eq.'con')then
            icnvel = 0
        else if(string(1:3).eq.'VAR' .or. string(1:3).eq.'var')then
            icnvel = 1
        endif
c-----
c       get lines 8 through 11
c-----
        do 900 i=8,11
            read(lun,'(a)')string
  900   continue
c-----
c       get model specifically for 1-D flat isotropic
c-----
c-----
c       get comment line
c-----
        read(lun,'(a)')string
        mmax = 0
        refdep = 0.0
        irefdp = 0
        if(iiso.eq.0)then
 1000       continue
            j = mmax +1
                read(lun,*,err=9000,end=9000)d(j),a(j),b(j),
     1              rho(j),qa(j),qb(j),etap(j),etas(j),
     2              frefp(j),frefs(j)
                if(d(j).lt.0.0)then
                    d(j) = -d(j)
                    refdep = refdep + d(j)
                    irefdp = j
                endif
            mmax = j
            go to 1000
 9000       continue
        endif
    1   format(' LAYER             H      P-VEL     S-VEL   DENSITY  ')
    2   format(' ',i5,5x,4f10.3)
    3   format(' ','-SURFACE ','- - - - - ','- - - - - ',
     1      '- - - - - ','- - - - - -')
        if(mmax.gt.0)then
            if(listmd)then
            ierr = 0
            write(LOT,1)
            do 2000 i=1,mmax
                write(LOT,2)
     1              i,d(i),a(i),b(i),rho(i)
                if(i.eq.irefdp)write(LOT,3)
 2000       continue
            endif
        else 
            ierr = -3
            write(LER,*)'Error in model file'
        endif
        if(lun.ne.LIN)close (lun)
        return
        end
        function lgstr(str)
c-----
c       function to find the length of a string
c       this will only be used with file system path names
c       thus the first blank 
c       indicates the end of the string
c-----
        implicit none
        character*(*) str
        integer*4 lgstr
        integer n, i
        n = len(str)
        lgstr = 1
        do 1000 i=n,1,-1
            lgstr = i
            if(str(i:i).ne.' ')goto 100
 1000   continue
  100   continue
        return
        end
