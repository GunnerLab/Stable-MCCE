        subroutine phierg(nit1,natom,rtemp,cutoff,nqass)
c analyse phimap energies etc
c
c-------------------------------------------------
        include "qdiffpar4.h"
        include "qlog.h"
c-------------------------------------------------
c
        
      real phimap(igrid,igrid,igrid),rad3(natom)
      real*4 qmap(igrid,igrid,igrid)
      real*4 debmap(igrid,igrid,igrid)
      integer iepsmp(igrid,igrid,igrid,3)
      integer nepsmp(igrid,igrid,igrid,3)
        
        real atmcrg(4,nqass)
        
        dimension nclass(7,2,2),indx(3)
        dimension qass(4,natmax),inear(6)
        real*4 cputime,start,fin
        real*8 terg,qfix,trhophi,tosm,tmp
        character*8 hour
        logical donon
        
c      data twopi / 6.2832 /
c      data fpi / 12.566 /
c      data fact  / 6.0e-4  /   !converts moles to ions/cubic 
c        angstrom
c                                ! ==Na,l/moles,Ang*3
c        data conv / 1.78e-3  /         !converts integration to units of e
c                                !==erg,cm,e,e/kT,Ang,esu,esu
c-------------------------------------------------------
        start = cputime(0.0)
        call time(hour)
        print *,' '
        write(6,*)'analysing phimap at: ',hour
c-------------------------------------------------------


        do k=1,igrid
          do j=1,igrid
            do i=1,igrid              
              do m=1,3
                if(iepsmp(i,j,k,m).gt.natom) then
                  nepsmp(i,j,k,m)=1
                else
                  nepsmp(i,j,k,m)=0
                end if                
              end do
            end do
          end do
        end do
        
        
       pi = acos(-1.0)
       twopi = 2.0 * pi
       fpi = 4.0 * pi
       fact = 6.022136736e-4

        midg = (igrid+1)/2
c
c total of q.phi, includes grid energy
c
c        pi = 355./113.
        donon = .false.
        if ( nit1.ge.0 .and. rionst.gt.0.0 ) donon=.true. 

c        nqmap=0
c        nphimap=0
        
        terg = 0.0
        qfix = 0.0
        do k = 1,igrid
          do j = 1,igrid
            do i = 1,igrid
                qfix = qfix + qmap(i,j,k)
                terg = terg + qmap(i,j,k) * phimap(i,j,k)
c                if(qmap(i,j,k).eq.0.) nqmap=nqmap+1
c                if(phimap(i,j,k).eq.0.) nphimap=nphimap+1
            end do
          end do
        end do

c        print *,'number of times when qmap is 0 ',nqmap
c        print *,'number of times when phimap is 0 ',nphimap

        terg = terg / 8.0 / pi / scale
        print *,'0.5*sum(phi.q) (kT) =  ',terg

        qfix = qfix / 4.0 / pi / scale
        print *,'For net fixed charge of : ',qfix
c
c integrate normal field at edge of box to get net charge 
c
        eint = 0.d0
        darea = 1.0 / scale**2
        do i = 2,igrid-1
          do k = 2,igrid-1
            eint = eint + (phimap(k,1,i) - phimap(k,3,i))
            eint = eint + (phimap(k,igrid,i) - phimap(k,igrid-2,i))
            eint = eint + (phimap(1,k,i) - phimap(3,k,i))
            eint = eint + (phimap(igrid,k,i) - phimap(igrid-2,k,i))
            eint = eint + (phimap(k,i,1) - phimap(k,i,3))
            eint = eint + (phimap(k,i,igrid) - phimap(k,i,igrid-2))
          end do
        end do
        qint = - eint * epsout * darea * scale / fpi / 2.0
        print *,'net charge from gauss integral: ',qint
        if (rionst.gt.0.0) then
c
c initialize energy density integrals
c

        vol = 1.0 / scale**3
c cal z range indicies
          zlw = (2-midg) / scale + oldmid(3)
          zup = (igrid-1-midg) / scale + oldmid(3)
          print *,'integrating over zrange: ',zlw,zup

          k1=0
          
          tions = 0.0
          trhophi = 0.0
          tedotd = 0.0
          tosm = 0.0
          tcpos = 0.0
          tcneg = 0.0
*          do k = 2,igrid-1
*            do j = 2,igrid-1
*              do i = 2,igrid-1
          do k = 1,igrid
            do j = 1,igrid
              do i = 1,igrid
                if (debmap(i,j,k).gt.0.0) then
                  pot = phimap(i,j,k)
                  k1=k1+1
c charge concentrations 
                    if(donon)then
                      cpos = exp(-max(pot,-cutoff)) - 1.0
                      cneg = exp(min(pot,cutoff)) - 1.0
                    else
                      cpos = -pot 
                      cneg = pot
                    end if
                    cnet = cpos - cneg
                    tcpos = tcpos + cpos
                    tcneg = tcneg + cneg
                    tions = tions + cpos - cneg
c entropy contributions
                    trhophi = trhophi + pot*cnet
c osmotic pressure term
                    tosm = tosm + cpos + cneg
                end if
            end do   
          end do      
          end do
                

          do k = 2,igrid-1
            do j = 2,igrid-1
              do i = 2,igrid-1                
c
c electrostatic stress term (E.D)
c
                iin = (nepsmp(i,j,k,1) + nepsmp(i,j,k,2) 
     &               + nepsmp(i,j,k,3) + nepsmp(i-1,j,k,1)
     &               + nepsmp(i,j-1,k,2) + nepsmp(i,j,k-1,3))
c use eps map
                if(iin.eq.0) then
                    fldx = (phimap(i+1,j,k) - phimap(i-1,j,k)) / 2.0
                    fldy = (phimap(i,j+1,k) - phimap(i,j-1,k)) / 2.0
                    fldz = (phimap(i,j,k+1) - phimap(i,j,k-1)) / 2.0
                    tedotd = tedotd + (fldx**2 + fldy**2 + fldz**2)
                end if
            end do   
          end do      
          end do
          
          write(6,*) "Points have salt ",k1

          tions = tions * rionst * vol * fact
          trhophi = trhophi * rionst * fact * vol
          tedotd = tedotd * vol * scale * scale * epsout / fpi
          tosm = tosm * rionst * vol * fact
          tcpos = tcpos * rionst * vol * fact
          tcneg = tcneg * rionst * vol * fact
          ergion = tedotd / 2.0 - trhophi - tosm
          qnrat = tcneg / abs(tions)
          qprat = tcpos / abs(tions)
  
          calpkt = 1.98720 * rtemp / 1000.0
          print *,'kT = ',calpkt,' kcal at ',rtemp,'K'
          print *,' '
          print *,'Q net                       : ',tions
          print *,'Q+ excess,ratio             : ',tcpos,qprat
          print *,'Q- excess,ratio             : ',tcneg,qnrat
          print *,'rho.phi or -TdS          (kT): ',trhophi
          print *,'E.D                     (kT): ',tedotd
          print *,'dPI                     (kT): ',tosm
          print *,'E.D/2 - rho.phi - dPI   (kT): ',ergion
          print *,' '

c The salt part of eq. (24) in the Sharp & Honig paper:

        write(6,*) ' '
        write(6,*) ' - 0.5 * rho.phi - dPI =',- 0.5 * trhophi - tosm
        write(6,*) ' '
        fin = cputime(start)
        print *,'Salt analysis took ',fin,' sec'
        start = cputime(0.0)
        end if

        tmp = terg - 0.5*trhophi-tosm
        write(6,*) ' '   
        write(6,*) 
     & 'Calculate the energy using a subroutine from QNIFFT'
        write(6,*) ' linear_energy = ',terg
        write(6,*) ' non_linear_energy = ',tmp*calpkt,' (kcal/mol)'
        write(6,*) ' '                
                
c
c collect assigned charge
c
        nq = 0
        qfix = 0.0
        rmidg = (igrid+1)/2.
        do i = 1,natom
          if(atmcrg(4,i).ne.0.0)then
            nq = nq + 1
            do k = 1,3
                qass(k,nq) = atmcrg(k,i)
c                qass(k,nq) = (atmcrg(k,i)-rmidg)/scale+oldmid(k)
            end do
            qfix = qfix + atmcrg(4,i)
            qass(4,nq) = atmcrg(4,i)
c          print *,qass(4,nq),qass(1,nq),qass(2,nq),
c            qass(3,nq)
          end if
        end do
        print *,'total atomic charge: ',qfix
c
c coulombic energy
c
        ergc = 0.0
        do i = 1,nq
          do j = i+1,nq
            r2 = (qass(1,i)-qass(1,j))**2 + (qass(2,i)-qass(2,j))**2
     &    + (qass(3,i)-qass(3,j))**2
c            print *,'r2: ',r2
            if(r2.gt.1.e-6)ergc = ergc + qass(4,i)*qass(4,j)/sqrt(r2)
          end do
        end do
c        print *,'epsin: ',epsin
        ergc = ergc/epsin*scale
        print *,'Coulombic energy in uniform Epsin (kT): ',ergc
        fin = cputime(start)
        print *,'Coulombic energy took ',fin,' sec'
        start = cputime(0.0)
c
c surface charge reaction field energy
c
        ergs = 0.0
        qsurf = 0.0
        sfact = 0.25 / pi / scale / epkt
        nmove = 0
        nsurf=0
        nsurfq = 0
        do k = 2,igrid-1
          do j = 2,igrid-1
            do i = 2,igrid-1
              iin = (nepsmp(i,j,k,1) + nepsmp(i,j,k,2) 
     &             + nepsmp(i,j,k,3) + nepsmp(i-1,j,k,1)
     &             + nepsmp(i,j-1,k,2) + nepsmp(i,j,k-1,3))
              if ( iin.ne.0 .and. iin.ne.6 ) then
                  if (qmap(i,j,k).ne.0.0) then
                    nsurfq = nsurfq+1
                  else
                    nsurf = nsurf+1
                  end if
                  srfq = 6.0 * phimap(i,j,k) -
     &     (phimap(i+1,j,k) + phimap(i-1,j,k)
     &          + phimap(i,j+1,k) + phimap(i,j-1,k)
     &          + phimap(i,j,k+1) + phimap(i,j,k-1) + qmap(i,j,k) 
     &          / epsin)
                   qsurf = qsurf + srfq
c
c correction position of surface charge using position of nearest 
c atom
c ixepsmp2 contains indices of atoms forming surface
c
*                  inear(1) = ixepsmp2(i,j,k,1)
*                  inear(2) = ixepsmp2(i,j,k,2)
*                  inear(3) = ixepsmp2(i,j,k,3)
*                  inear(4) = ixepsmp2(i-1,j,k,1)
*                  inear(5) = ixepsmp2(i,j-1,k,2)
*                  inear(6) = ixepsmp2(i,j,k-1,3)
                  iat = 0
                  dmin2 = 1.e6
                  srfx = i
                  srfy = j
                  srfz = k
                  do i1 = 1,6
*                    if (inear(i1).ge.1 .and. inear(i1).le.natom) then
*                      dist2 = (atmcrg(1,inear(i1))-srfx)**2 
*     &                      + (atmcrg(2,inear(i1))-srfy)**2 
*     &                      + (atmcrg(3,inear(i1))-srfz)**2
*                      if (dist2.lt.dmin2) then
*                          dmin2 = dist2
*                          iat = inear(i1)
*                        end if
*                    end if
                  end do
c                  print *,'srfxyz: ',srfx,srfy,srfz
                  if (iat.ne.0) then     !found a neighbor atom
                    if (dmin2.gt.1.e-4) then
                      nmove = nmove + 1
c srat = atrad(iat) * scale / sqrt(dmin2) assuming atrat is rad3!!!
                   srat = rad3(iat) * scale / sqrt(dmin2)
c                  print *,'iat,rad,dmin2,srat: ',iat,atrad(iat),
c                   dmin2,srat
c                  print *,'atxyz: ',atmcrg(1,iat),atmcrg(2,iat),
c                   atmcrg(3,iat)
                      srfx = atmcrg(1,iat) + srat * (srfx-atmcrg(1,iat))
                      srfy = atmcrg(2,iat) + srat * (srfy-atmcrg(2,iat))
                      srfz = atmcrg(3,iat) + srat * (srfz-atmcrg(3,iat))
                    end if
                  end if
c                  print *,'srfxyz after correction: ',srfx,srfy,srfz
c
c loop thru charges and calc reaction energy
c
                  do i1 = 1,nq
                        r2 = (qass(1,i1)-srfx)**2 + (qass(2,i1)-srfy)**2 
     &             + (qass(3,i1)-srfz)**2
                  if (r2.gt.1.d-6) ergs = ergs + qass(4,i1)*
     &              srfq/sqrt(r2)
                  end do
              end if
          end do   
        end do      
        end do
        qsurf = qsurf*sfact
        qsurfe = (epsin-epsout)*qfix/epsin/epsout/epkt
        print *,'surface points readjusted: ',nmove
        print *,'surface points with, without fixed charge: ',
     &    nsurfq,nsurf
        print *,'expected, actual induced surface charge: ',qsurfe,
     &    qsurf
        ergs = ergs*epkt*sfact*scale/2.d0
        print *,'Dielectric rxn field energy from surface Q (kT): ',
     &    ergs
        fin = cputime(start)
        print *,'Dielectric Reaction energy took ',fin,' sec'

c The coulombic part is the sum of unifor dielectric and reaction 
c field
c energy. The salt part is given by eq. (24) in the Sharp & Honig 
c paper.

C        write(6,*) ' '        
C        write(6,*) ' linear_energy =',ergc + ergs
C        write(6,*) ' non_linear_energy =',ergc+ergs-0.5*trhophi-tosm
C        write(6,*) ' '                
                
        start = cputime(0.0)
        if(iconc)then
          print *,'converting potentials to Molar concentrations'
          do k = 1,igrid
            do j = 1,igrid
              do i = 1,igrid
                if(debmap(i,j,k).gt.0.)then
                    cnet = -rionst * 2.d0 * sinh(phimap(i,j,k))
                  else
                    cnet = 0.0
                end if
                  phimap(i,j,k) = cnet
              end do
            end do
          end do
        end if
        
        return
        end
