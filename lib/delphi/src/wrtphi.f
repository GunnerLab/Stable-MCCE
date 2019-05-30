	subroutine wrtphi
c
cNB uses phimap3 as a temp storage for phimap in case
c want to write not potentials bu salt concentraions
c
	include "qdiffpar4.h"
	include "qlog.h"
        parameter (mgrid=65)
c
      real phimap(igrid,igrid,igrid)
	character*80 filnam
	character*10 nxtlbl
      character*20 uplbl
      character*16 botlbl
c b+++++++++++++++++++++++++++for back compatibility with single precision phi readers
      real*4 phimap4(1),phimap3(ngp)
	real*4 scalesingle,scalesingle1,oldmidsingle(3),oldmidsingle1(3)
	real*4 minim,maxim,somma,average
      integer igrid1      

	scalesingle=scale
      oldmidsingle(1)=oldmid(1)
      oldmidsingle(2)=oldmid(2)
      oldmidsingle(3)=oldmid(3)


      if (realsiz.ne.4.and.phifrm.ne.2) then
	   i_phimap4=memalloc(i_phimap4,realsiz,igrid*igrid*igrid)
	   i=1
         do iz=1,igrid
         do iy=1,igrid
         do ix=1,igrid
            phimap4(i)=phimap(ix,iy,iz)
	      i=i+1
         end do
         end do
         end do
      end if
	  write(6,*) phimap(5,5,5)
	  write(6,*) phimap4(5+(5-1)*igrid+(5-1)*igrid*igrid)

c e++++++++++++++++++++++++++++
c
	if(iconc) then
        i_phimap3= memalloc(i_phimap3,realsiz,ngp)
	i=1
	do iz=1,igrid
	do iy=1,igrid
	do ix=1,igrid
	phimap3(i)=phimap(ix,iy,iz)
	i=i+1
	end do
	end do
	end do
	call phicon
	end if
c
	if(ibios) then
c
c write phimap in insight format
c

	  open(14,file=phinam(:philen),form="unformatted")
	  filnam = ' '
	  inquire(14,name = filnam)
	  write(6,*)'potential map written in INSIGHT format to file'
	  write(6,*)filnam
	  write(6,*)'  '
	  ivary = 0
	  nbyte = 4
	  intdat = 0
	  xang = 90.
	  yang = 90.
	  zang = 90.
	  intx = igrid - 1
	  inty = igrid - 1
	  intz = igrid - 1
	  xmax = 0.
	  do 9040 k = 1,3
	    temp = abs(oldmid(k))
	    xmax = max(xmax,temp)
9040	  continue
	  range = (igrid-1.)/(2.*scale)
	  extent = range + xmax
	  xstart = (oldmid(1)-range)/extent
	  ystart = (oldmid(2)-range)/extent
	  zstart = (oldmid(3)-range)/extent
	  xend = (oldmid(1)+range)/extent
	  yend = (oldmid(2)+range)/extent
	  zend = (oldmid(3)+range)/extent
        write(14)toplbl
	  write(14)ivary,nbyte,intdat,extent,extent,extent,
     1  xang,yang,zang,xstart,xend,ystart,yend,zstart,
     1  zend,intx,inty,intz
         print*,' in wrtphi, realsiz,toplbl=', realsiz,toplbl
	  if (realsiz.ne.4) then
          call wrtphimap(igrid,phimap4,1)
	  else
	    do 9044 k = 1,igrid
	    do 9043 j = 1,igrid
		  write(14)(phimap(i,j,k),i=1,igrid)
9043	    continue
9044	    continue
	  end if
 
c
	elseif(phifrm.eq.2)then
c GRASP phimap - output a 65^3 grid and leave out ngrid spec.
c
          write(6,*)'  '
          write(6,*)'writing potential map in GRASP format'
          write(6,*)'  '
c
          i_phimap4=memalloc(i_phimap4,realsiz,mgrid*mgrid*mgrid)
          call expand(mgrid)
          open(14,file=phinam(:philen),form="unformatted")
          filnam = ' '
          inquire(14,name = filnam)
          write(6,*)'potential map written to file'
          write(6,*)filnam
          write(6,*)'  '
          if(iconc.and.(rionst.ne.0)) then
          nxtlbl="concentrat"
          else
          nxtlbl="potential "
          end if
        write(14)'now starting phimap '
        write(14)nxtlbl,toplbl
c        write(14)phimap4
        call wrtphimap(mgrid,phimap4,0)
        write(14)' end of phimap  '
        write(14)scalesingle,oldmidsingle
        close(14)
          i_phimap4=memalloc(i_phimap4,0,0)
c b++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	elseif(phifrm.eq.3)then
c CCP4 phimap - output 
c
	  write(6,*)'  '
	  write(6,*)'writing potential map in CCP4 format'
	  write(6,*)'  ' 
c       '(A20)'
	  open(14,file=phinam(:philen),form="unformatted")
	    filnam = ' '
	    inquire(14,name = filnam)
	    write(6,*)'potential map written to file'
	    write(6,*)filnam
  	    write(6,*)'  '
c NC,NR,NS,MODE,NCSTART,NRSTART,NSSTART,NX,NY,NZ
          write(14)igrid,igrid,igrid,2,1,1,1,igrid-1,igrid-1,igrid-1
c X Y Z lengths, alpha
	    write(14)(igrid-1)/scale,(igrid-1)/scale,(igrid-1)/scale,90.0
c beta, gamma,MAPC,MAPR,MAPS 
		write(14)90.0,90.0,1,2,3
c 
		minim=10000.0
	    maxim=-10000.0
	    if (realsiz.ne.4) then
	        i=1
			do iz=1,igrid
			do iy=1,igrid
			do ix=1,igrid
			  if (phimap4(i).lt.minim) minim=phimap4(i)
			  if (phimap4(i).gt.maxim) maxim=phimap4(i)
			  somma=somma+phimap4(i)
	          i=i+1
			end do
			end do
			end do
	        average=somma/(igrid*igrid*igrid)
	    else
			do iz=1,igrid
			do iy=1,igrid
			do ix=1,igrid
				if (phimap(ix,iy,iz).lt.minim) minim=phimap(ix,iy,iz)
				if (phimap(ix,iy,iz).gt.maxim) maxim=phimap(ix,iy,iz)
				somma=somma+phimap(ix,iy,iz)
			end do
			end do
			end do
              average=somma/(igrid*igrid*igrid)
	    end if
c amin, amax, amean,ispg,nsymbt,LSKFLG
		write(14)minim,maxim,average,1,1,0
c skwmat(3,3)
          write(14)0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
c skwtrn(3),future use
          write(14)0.0,0.0,0.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
          write(14)'MAP '
c MACHST,ARMS,NLABL
          write(14)1,0.0,0
	    do ii=1,200
	      write(14)0
	    end do
	    if (realsiz.ne.4) then
	      call wrtphimap(igrid,phimap4,0)
	    else
             write(14)phimap
	    end if
          close(14)

	elseif(phifrm.eq.4)then
c diffrential phimap - 
c
          write(6,*)'  '
          write(6,*)'reading previous.phi file'
	    write(6,*) realsiz, phimap(5,5,5)
	    write(6,*) phimap4(5+(5-1)*igrid+(5-1)*igrid*igrid)
          write(6,*)'  '
c
          open(14,file="previous.phi",form="unformatted")
             read(14)uplbl
	       read(14)nxtlbl,toplbl
             write(6,*)uplbl,nxtlbl,toplbl
             call rdphimap(igrid,phimap4)
	       read(14)botlbl
		   read(14)scalesingle1,oldmidsingle1,igrid1
	       write(6,*)botlbl,scalesingle1,oldmidsingle1
          close(14)
          if ((scalesingle1.ne.scalesingle).or.(oldmidsingle1(1).ne.
     &      oldmidsingle(1)).or.(oldmidsingle1(2).ne.oldmidsingle(2)).
     &      or.(oldmidsingle1(3).ne.oldmidsingle(3)).or.(igrid1.ne.
     &      igrid)) then
	        write(6,*) 'Error: the two potential maps do not mach'
	        return
	    endif
          write(6,*)phimap4(5+(5-1)*igrid+(5-1)*igrid*igrid)
	    if (realsiz.ne.4) then
             call submaps(igrid,phimap,phimap4,2)
	    else
             call submaps(igrid,phimap,phimap4,1)
	    end if

	    open(14,file=phinam(:philen),form="unformatted")
	    filnam = ' '
	    inquire(14,name = filnam)
	    write(6,*)'potential map written to file'
	    write(6,*)filnam
  	    write(6,*)'  '
  	    if(iconc.and.(rionst.ne.0)) then
  	      nxtlbl="concentrat"
  	    else
	      nxtlbl="potential "
	    end if
c         the following is uplabel, 20 char
          write(14)'now starting phimap '
          write(14)nxtlbl,toplbl
	    if (realsiz.ne.4) then
             call wrtphimap(igrid,phimap4,0)
	    else
             write(14)phimap
	    end if
c         this is botlbl
          write(14)' end of phimap  '
          write(14)scalesingle,oldmidsingle,igrid
          close(14)

c e++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      else
	  write(6,*)'  '
	  write(6,*)'writing potential map in DELPHI format'
	  write(6,*) realsiz, phimap(5,5,5)
	  write(6,*) phimap4(5+(5-1)*igrid+(5-1)*igrid*igrid)
	  write(6,*)'  '
c
	  open(14,file=phinam(:philen), form='formatted')
	    filnam = ' '
	    inquire(14,name = filnam)
	    write(6,*)'potential map written to file'
	    write(6,*)filnam
  	    write(6,*)'  '
  	    if(iconc.and.(rionst.ne.0)) then
  	      nxtlbl="concentrat"
  	    else
	      nxtlbl="potential "
	    end if
c         the following is uplabel, 20 char
          write(14,'(a20)')'now starting phimap '
          write(14,'(a10,a60)')nxtlbl,toplbl
             write(14,'(7f15.7)')phimap
c         this is botlbl
          write(14,'(a16)')' end of phimap  '
          write(14,'(f15.7,3f15.7,i10)')scalesingle,oldmidsingle,igrid
          close(14)
	end if
	
      if(realsiz.ne.4.and.phifrm.ne.2) i_phimap4=memalloc(i_phimap4,0,0)
c
c b ++++++debug++++++++++++
c       open(52,FILE='phiwalt')
c       do ix=1,65
c            do iz=1,65
c              write(52,*)phimap(ix,33,iz)
c            end do
c       end do
c       close (52)
c  e +++++++++++++++++++++++

	if(iconc) then
	i=1
	do iz=1,igrid
	do iy=1,igrid
	do ix=1,igrid
	phimap(ix,iy,iz)=phimap3(i)
	i=i+1
	end do
	end do
	end do
        i_phimap3= memalloc(i_phimap3,0,0)
	end if
c
	return
	end

      subroutine wrtphimap(mgrid,phimap,opt)
      real*4 phimap(mgrid,mgrid,mgrid)
	integer opt

	 if (opt.eq.1) then
	   do 9041 k = 1,mgrid
	   do 9042 j = 1,mgrid
	     write(14)(phimap(i,j,k),i=1,mgrid)
9042	   continue
9041	   continue
       else
         write(14)phimap
	 endif
       return
      end

      subroutine rdphimap(mgrid,phimap)
      real phimap(mgrid,mgrid,mgrid)
         read(14)phimap
	   write(6,*) 'inside rdphimap',phimap(5,5,5)
       return
      end

      subroutine submaps(mgrid,phimap,phimap4,opt)
      integer mgrid,opt
	real phimap(mgrid,mgrid,mgrid)
      real phimap4(mgrid,mgrid,mgrid)

      if(opt.eq.1) then
	  do iz=1,mgrid
	  do iy=1,mgrid
	  do ix=1,mgrid
	    phimap(ix,iy,iz)=phimap(ix,iy,iz)-phimap4(ix,iy,iz)
	  end do
	  end do
	  end do
	else
	  do iz=1,mgrid
	  do iy=1,mgrid
	  do ix=1,mgrid
	    phimap4(ix,iy,iz)=phimap(ix,iy,iz)-phimap4(ix,iy,iz)
	  end do
	  end do
	  end do
	endif

        return
      end

        subroutine wrtphiForGUI
c
        include "qdiffpar4.h"
        include "qlog.h"
c
        real phimap(igrid,igrid,igrid)
c b++++++++++++++walter++++write potential map for the GUI
c       open(52,FILE="phimap",form='unformatted')
c       write(52)phimap
c       close(52)
c e+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        return
        end


