	subroutine encalc_qnifft(icount1b,nqass,natom,ibnum,nmedia,nqgrd
     &    ,donon,cutoff) 
c 
	include "qdiffpar4.h" 
	include "qlog.h" 
c 
	real phimap(igrid,igrid,igrid) 
	real debmap(igrid,igrid,igrid)
	logical*1 idebmap(igrid,igrid,igrid)
      logical donon
	dimension gchrg(icount1b) 
	dimension rad3(natom),xn2(3,natom),scspos(3,ibnum) 
	integer gchrgp(3,icount1b),ibgrd(3,ibnum),nqgrd 
	logical ido 
c b++++++++++++++ 
        integer igridout 
        real ergs,ergas,ergc         
        real*8 ergg,ergnl,calpkt,trhophi,tosm                
c 
	if(inrgwrt) open(42,file=nrgnam(:nrgfrm)) 
c 
	if(loga)  then 
	call anagrd(icount1b,epsin*epkt,erga,scale) 
	write(6,*) 'analytic grid energy is        ',erga,' kt' 
	if(inrgwrt) write(42,*) 'analytic grid energy is',erga,' kt' 
	end if 
c 
	if(logg) then 
	  ergg=0.0 
          limx1=2+bufz(1,1) 
	  limx2=igrid-1-bufz(2,1) 
	  limy1=2+bufz(1,2) 
	  limy2=igrid-1-bufz(2,2) 
	  limz1=2+bufz(1,3) 
	  limz2=igrid-1-bufz(2,3) 
	  do 587 i=1,icount1b 
	    ix=gchrgp(1,i) 
	    iy=gchrgp(2,i) 
	    iz=gchrgp(3,i) 
	    ido=.true. 
	    if((ix.lt.limx1).or.(ix.gt.limx2)) ido=.false. 
	    if((iy.lt.limy1).or.(iy.gt.limy2)) ido=.false. 
	    if((iz.lt.limz1).or.(iz.gt.limz2)) ido=.false. 
	    if(ido) ergg=ergg + phimap(ix,iy,iz)*gchrg(i) 
587	  continue 
	  ergg=ergg/2.0 
	  
	  calpkt = 1.98720 * 298.15 / 1000.0
	  write(6,*) ' ' 
	  write(6,*) 'total grid energy          :     ',ergg,' kt' 
	  if(inrgwrt) write(42,*)'total grid energy: ',ergg,' kt' 
	end if 
c 
	if(logg.and.loga) then 
	  write(6,*) 'difference energy, in kt, is',(ergg-erga) 
	  write(6,*) 'difference energy, in kcals, is',(ergg-erga)*0.6 
	end if 
c 
c b+++++++++++++++++++++++++++++++++++++w Oct 2000 
	if(irea.or.logs.or.lognl.or.logas.or.isen.or.isch) then 
          ergs=0.0 
          ergas=0.0 
          ergnl=0.0 
          ergest=0.0 
c ergest=interaction energy of the solvent and the fixed charges 
	  if(diff) ibc=0 
	  iisitpot=0
	  if(isitpot) iisitpot=1
	  call react(nqass,icount1b,ibnum,ergs,ergas,natom,
     &      nmedia,nobject,iisitpot) 
	end if 
c 
	if(logc.and.(.not.logions.or..not.lognl)) then 
          ergc=0.0 
          if (logions) then 
c           linear case 
            ergest=0.0 
            call clbtot(nqass,ergest,ergc) 
            write(6,*)'solvent contribution to fixed charges' 
            write(6,*)'respectively inside and outside the cube :', 
     &ergest,'kt',ergestout,'kt' 
            write(6,*)'total ionic direct contribution :',ergest+ 
     &ergestout,'kt' 
          else 
            if (nmedia.eq.1) then 
 	      call clb(nqass,ergc) 
              ergc=ergc/epsin 
            else 
              call clbmedia(nqass,ergc) 
            end if 
          end if 
	  write(6,*) 'coulombic energy :              ',ergc,' kt' 
	  if(inrgwrt) write(42,*) 'total coulombic energy:',ergc,' kt' 
        end if 
 
*        if (lognl) then 

* Joe modified
* there are some problems with the subroutine nlener
* Thus we use the formula from QNIFFT

*          call nlener(ergnl,igridout,phimap,idebmap) 

          trhophi=0.0
          tosm=0.0          
          k1=0
          
          do k = 1,igrid
            do j = 1,igrid
              do i = 1,igrid
                if (debmap(i,j,k).gt.0.0) then
                  pot = phimap(i,j,k)                  
c charge concentrations 
                    if(donon)then
                      if((pot.gt.cutoff).or.(pot.lt.-cutoff)) then
                        k1=k1+1
                      endif
                      cpos = exp(-max(pot,-cutoff)) - 1.0
                      cneg = exp(min(pot,cutoff)) - 1.0
                    else
                      cpos = -pot 
                      cneg = pot
                    end if
                    cnet = cpos - cneg
c entropy contributions
                    trhophi = trhophi + pot*cnet
c osmotic pressure term
                    tosm = tosm + cpos + cneg
                end if
            end do   
          end do      
          end do
          
          vol=1.0/scale**3          
          fact=6.022136736e-4
          
          trhophi = trhophi * rionst * fact * vol
          tosm = tosm * rionst * vol * fact     

          write(6,*)'Number of abs(phi) values over cutoff (50)',k1
          
          write(6,*)'rho*phi/2 term in solution :      ',
     &      trhophi/2.0,'kt'  
          write(6,*)'osmotic pressure term      :      ',
     &      tosm,'kt' 
          
          ergnl = -trhophi/2.0-tosm
                    
          write(6,*) 
     &  'Calculate the energy using a subroutine from DELPHI'
	  write(6,*) 'linear_energy = ',ergg*calpkt,' (kcal/mol)' 
	  write(6,*) 'non_linear_energy = ',(ergg+ergnl)*calpkt,
     &      ' (kcal/mol)' 
 
          ergc=0.0 
          ergest=0.0 
          if (logions) then 
            call clbnonl(nqass,ergc,ergest,igridout) 
          write(6,*)'direct ionic contrib. inside the box:',ergest,' kt' 
          write(6,*) 'coulombic energy:                     ',ergc,' kt' 
          if(inrgwrt) write(42,*) 'total coulombic energy:',ergc,' kt' 
          end if 
*        end if 
 
        if (logs.and.logions) then 
          write(6,*)'Energy arising from solvent and boundary pol.', 
     &ergnl+ergs+ergest+ergestout,' kt' 
        end if 
 
        ergtot=ergnl+ergc+ergs+ergest+ergestout 
        write(6,*)'All energy terms but grid and self_react.:',ergtot, 
     &'kt' 
        if(inrgwrt) write(42,*)'total energy (everything calculated  
     & but grid and self_reaction energies: ',ergtot,'kt' 
c e+++++++++++++++++++++++++++++++++++++ 
c 
	if(logas.and.loga.and.logg) then 
	write(6,*) "excess grid energy= ",ergg-ergas-erga 
	end if 
c 
	finish=cputime(start) 
	write(6,*) 'energy calculations done at',finish 
c 
	if(inrgwrt) close(42) 
c 
	return 
	end 
