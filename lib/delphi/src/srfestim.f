	subroutine srfestim(dx,dy,dz,R,scale)
c+++++++++++added by walter 1 Feb 2005
c definizione della funzione di area di una sfera intersecante un cubo.
      
	integer icont,sqrdim,identro,icontmax
      real dx,dy,dz,rho,rhoinv,deltateta,sindteta,fact,ca,cb,sa,sb
	real ll,passod,passoR,dmin,dmax,tetacen,phicen,deltaphi
	real alpha,beta,ii,x,y,z,rapp,area,acalotta,scale
	real pi,step,r1x,r1y,r1z,r2x,r2y,r3x,r3y,r3z,R
         
      pi=3.1416
	sqrdim=500
c     assumo lato cubo pari a 1
      ll=1./scale

c     rho= distanza tra centro cubo e centro sfera 
      rho=sqrt(dx*dx+dy*dy+dz*dz)
      rhoinv=1./rho
c     teta scandisce la verticale e phi il piano xy
c     mi devo muovere in un intorno del cubetto
c     calcolo i parametri angolari del centro del cubetto 
      tetacen=pi/2.
	phicen=0.
      deltaphi=1.3*ll*rhoinv
	deltateta=1.3*ll*rhoinv
      sindteta=sin(deltateta)
      fact=4*deltaphi*sindteta
c	 *sin(tetacen)
c      materiale per i cambiamenti di base

      if ((dx.eq.0.).and.(dy.eq.0.)) then
        alpha=0.
      else                
        alpha=-asin(dy/sqrt(dx*dx+dy*dy))
      end if
      beta=-asin(dz*rhoinv)
                        
      ca=cos(alpha) 
	sa=sin(alpha)
      cb=cos(beta)  
	sb=sin(beta)                

	r1x=cb*ca
	r1y=cb*sa
	r1z=sb
	r2x=-sa
	r2y=ca	
	r3x=-sb*ca
	r3y=-sb*sa
	r3z=cb      

c     genero una griglia di valori a caso di cos(teta) e phi
      step=2./sqrdim
	icont=1
	identro=0      
      do cteta=-sindteta,sindteta,step*sindteta
	   steta=sqrt(1.-cteta**2)
	   do phi=-deltaphi,deltaphi,step*deltaphi
	     icont=icont+1	
c     ora devo vedere se x,y e z dello stesso punto sono compresi nel cubetto
           x=steta*(cos(phi)*r1x+sin(phi)*r1y)+cteta*r1z
	     if (abs(2.*(x*R-dx)).lt.ll) then
             y=steta*(cos(phi)*r2x+sin(phi)*r2y)
	       if (abs(2.*(y*R-dy)).lt.ll) then
               z=steta*(cos(phi)*r3x+sin(phi)*r3y)+cteta*r3z
	         if (abs(2.*(z*R-dz)).lt.ll) then
	           identro=identro+1
	         end if
	       end if
	     end if                
         end do
	end do
      icontmax=icont-1

      rapp=float(identro)/icontmax
c      if(identro*1000.lt.icontmax) then
c        write(6,*) 'pochi per statistica',rapp
c      end if
c	 %*sin(tetacen)
      area=rapp*R*R*fact
	end
