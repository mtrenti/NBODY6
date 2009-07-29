c	specimen input file
c-0.00152 0.00607 1000.0  0.004299 5 1000 2.5
c
      PROGRAM aking6
c	This version computes the King model first.  It also attempts to be
c       faster than aking5 by improving the computation of N-body initial
c       conditions
c
c	Program to make an anisotropic King model (cf. Heggie & Ramamani,
c	MNRAS, 272, 317-22, 1995.)
c	Written (thrown together) by d.c.heggie@uk.ac.ed
c	This program comes without warranty
c	
c	Unit 6 shows, first, the course of the iteration in which
c	the analytic model is computed, second, the density
c	profiles along the three axes (in program units), and, third,
c	the limiting radii along the three axes, in both program
c	and user units
c	Unit 8 produces a list of the analytic solution
c 	Unit 9 produces masses, coordinates and velocity
c	components, in the user's units
c       Unit 7 produces surface densities (along the x and z axes,
c       as viewed along the y-axis) if N<=1
      implicit double precision (a-h,o-z)
      PARAMETER(NE=5,M=51,NB=3,NCI=NE,NCJ=NE-NB+1,NCK=M+1,NSI=NE,
     *     NSJ=2*NE+1,NYJ=NE,NYK=M,nstarm = 50000)
      COMMON X(M),H,MM,N,C2,ANORM
      COMMON /PATH/ KMAX,KOUNT,DXSAV,XP(200),YP(10,200)
      double precision mtot,mstar,j22,ke,kevir
      common /params/ w0,alpha1,alpha3,rho(m),rhodot(m),rlag,yking(2,m),
     &     den,rking
      DIMENSION SCALV(NE),INDEXV(NE),Y(NE,M),C(NCI,NCJ,NCK),S(NSI,NSJ)
	parameter (nmax=10)
	dimension ystart(nmax),coord(nstarm,3)
c      data scalv/4*0.001d0,1.d0/
      data scalv/5*1.d0/
	external derivs,rkqc
      densty(z) = -sqrt(z)*(z+1.5d0)+0.75*sqrt(pi)*exp(z)*erf(sqrt(z))
      pi = 4.d0*atan(1.d0)
      ITMAX=100
      CONV=5.E-6
      H=1.d0/(M-1)
      WRITE (*,*) 'ENTER alpha1, alpha3, mtot, G, iseed, N, W0 in a',
     &     ' consistent set of units'	
c	The galactic tidal potential (with the convention that
c	the potential is negative deep inside a bound system) is
c	0.5(alpha1*x**2 + alpha3*z**2)
      READ (*,*)  alpha1, alpha3, mtot, G, iseed, N, W0
	if (alpha1.ge.0.d0) then
		write (6,*) 'alpha1 should be negative'
		stop
	endif
	nread = 0
	itmax = 100
	rmin = 0.2
	if (w0.gt.12.d0) then
		write (6,*) 'w0 too large'
		stop
	endif
	if (w0.lt.0.2d0) then
		write (6,*)'w0 too small'
		stop
	endif
        slowc = 1.d0
      do i = 1,5
         indexv(i) = i
      enddo
      ANORM=1.
      h = (w0-rmin)**2/(m-1)
      den = densty(w0)
	dxsav = 0.01d0*w0
	kmax = 200
	nvar = 2
        xstart = 0.02
c        xstart = 0.001
	ystart(1) = 2.d0*xstart/3.d0
	ystart(2) = -2.d0/3.d0
c
c       In the computation of the King model, y1 = r^2, y2 = d(y1)/dx
	x1 = w0 - xstart
	x2 = 0.d0
	tol = 1.d-6
	h1 = 0.01d0
	hmin = 0.0000000001d0
      call ODEINT(YSTART,NVAR,X1,X2,tol,H1,HMIN,NOK,NBAD,DERIVS,RKQC)
      DO 12 K=1,M
         X(K) = w0-xstart-sqrt(h*(k-1))
         if (x(k).ge.0.d0) then
            rho(k) = densty(x(k))/den
            rhodot(k) = rho(k) + x(k)**1.5d0/den
         else
            rho(k) = 0.d0
            rhodot(k) = 0.d0
         endif
	if (x(k).gt.xp(1)) then
		yking(1,k) = (w0-x(k))*yp(1,1)/(w0 - xp(1))
		yking(2,k) = yp(2,1)
	else
		i = 1
115		continue
		if ((x(k)-xp(i))*(x(k)-xp(i+1)).le.0.d0) then
			yking(1,k) = yp(1,i) + (yp(1,i+1)-yp(1,i))*
     &				(x(k)-xp(i))/(xp(i+1)-xp(i))
			yking(2,k) = yp(2,i) + (yp(2,i+1)-yp(2,i))*
     &				(x(k)-xp(i))/(xp(i+1)-xp(i))
		else
			i = i+1
			if (i.le.kount-1) then
				goto 115
			else
				write (6,*) 'failing interpolation',
     &					' of King values?',xp(kount),
     &                                  x(k)
                                yking(1,k) = yp(1,kount)
                                yking(2,k) = yp(2,kount)
			endif
		endif
	endif
         y(1,k) = -0.000004*(1-x(k)/w0)
         y(2,k) = 0.000004/w0
         y(3,k) = 0.000001*(1-x(k)/w0)
         y(4,k) = -0.000001/w0
         y(5,k) = log(0.000002d0)
   12 CONTINUE
      if (nread.eq.1) then
         do j = 1,m
            read (7,100) a,(y(i,j),i=1,5)
         enddo
      endif
      CALL SOLVDE(ITMAX,CONV,SLOWC,SCALV,INDEXV,NE,NB,M,Y,NYJ,NYK,
     *     C,NCI,NCJ,NCK,S,NSI,NSJ)
      do j = 1,m
         write (8,100) x(j),(yking(i,j),i=1,2),(y(i,j),i=1,5)
      enddo
  100 format (0pf10.5,1p7e10.2)
c	Now compute density profiles
      s0 = 1.5d0*(1.d0+alpha3/alpha1)
      s2x = -0.5d0*alpha3/alpha1 + 1.d0
      s2y = s2x - 1.5d0
      s2z = alpha3/alpha1 - 0.5d0
      eps = - exp(y(5,m))
      wyold = -1.d0
      wzold = -1.d0
	write (6,*) 'radius    density on x      on y      on z'
        fmax = 0.d0
      do i = 1,m
         r = sqrt(yking(1,i))
         wx = x(i) + y(1,i)*s0 + y(3,i)*s2x - eps*r**2*
     &        (2.d0/3.d0)*(s0/3.d0 + s2x)
         wy = x(i) + y(1,i)*s0 + y(3,i)*s2y - eps*r**2*
     &        (2.d0/3.d0)*(s0/3.d0 + s2y)
         wz = x(i) + y(1,i)*s0 + y(3,i)*s2z - eps*r**2*
     &        (2.d0/3.d0)*(s0/3.d0 + s2z)
         if (wx.ge.0.d0) then 
            rhox = densty(wx)/den
         else
            rhox = 0.d0
         endif
         if (wy.ge.0.d0) then 
            rhoy = densty(wy)/den
         else
            rhoy = 0.d0
         endif
         if (wz.ge.0.d0) then 
            rhoz = densty(wz)/den
         else
            rhoz = 0.d0
         endif
         write (6,101) r,rhox,rhoy,rhoz,rho(i)
         if (fmax.lt.r**2*rhox) fmax = r**2*rhox
  101    format (0pf10.5,1p4d15.5)
         if (wyold.gt.0.d0) ymax = (r*wyold - wy*
     &        sqrt(yking(1,i-1)))/(wyold-wy)
         if (wzold.gt.0.d0) zmax = (r*wzold - wz*
     &        sqrt(yking(1,i-1)))/(wzold-wz)
         wyold = wy
         wzold = wz
      enddo
      r = r + 0.1
   10 continue
      wc1 = -2.d0*yking(1,m)**1.5d0/yking(2,m)
      wc2 = -wc1/sqrt(yking(1,m))
      w0c1 = wc1*y(2,m)
      w0c2 = y(1,m) - w0c1/sqrt(yking(1,m))
      w2c1 = yking(1,m)**1.5d0*y(3,m)
      wx = wc1/r + wc2
     &     +s0*(w0c1/r + w0c2) + s2x*w2c1/r**3 - eps*r**2*
     &     (2.d0/3.d0)*(s0/3.d0 + s2x)
      wy = wc1/r + wc2
     &     +s0*(w0c1/r + w0c2) + s2y*w2c1/r**3 - eps*r**2*
     &     (2.d0/3.d0)*(s0/3.d0 + s2y)
      wz = wc1/r + wc2
     &     +s0*(w0c1/r + w0c2) + s2z*w2c1/r**3 - eps*r**2*
     &     (2.d0/3.d0)*(s0/3.d0 + s2z)
      if (wx.ge.0.d0) then 
         rhox = densty(wx)/den
      else
         rhox = 0.d0
      endif
      if (wy.ge.0.d0) then 
         rhoy = densty(wy)/den
      else
         rhoy = 0.d0
      endif
      if (wz.ge.0.d0) then 
         rhoz = densty(wz)/den
      else
         rhoz = 0.d0
      endif
      if (wyold.gt.0.d0) ymax = (r*wyold - wy*(r-0.1d0))/
     &     (wyold-wy)
      if (wzold.gt.0.d0) zmax = (r*wzold - wz*(r-0.1d0))/
     &     (wzold-wz)
      wyold = wy
      wzold = wz
      rhok = 0.d0
      write (6,101) r,rhox,rhoy,rhoz,rhok
         if (fmax.lt.r**2*rhox) fmax = r**2*rhox
      r = r + 0.1
      if (r.lt.rlag) goto 10
      fmax = 1.1d0
c      write (6,*) 'mass of King model, and anisotropic model ',wc1,
c     &     wc1+s0*w0c1
      write (6,*) 'radii (program units):',rlag,ymax,zmax
      rking = sqrt(yking(1,m))
      write (6,*) 'King radius (program units):',rking
      write (6,'(''Densities: King '',f10.5, '' Aniso '', f10.5,
     & '' W0 '',f10.3)') wc1*3/(4*pi*rking**3), 
     &     (wc1 + s0*w0c1)*3/(4*pi*rlag**3), W0
c
c      Next section computes surface density viewed parallel to y-axis
      if (n.le.1) then
         zero = 0.d0
         sigma0 = sigma(zero,zero,y)
         write (6,*) 'Computing surface densities...'
         do i = 1,250
            x1 = 10.d0**(0.01d0*(i-1))
            if (x1.le.rlag) then
               sigmax = sigma(x1,zero,y)/sigma0
               sigmaz = sigma(zero,x1,y)/sigma0
               write (7,*) x1,sigmax,sigmaz
            endif
         enddo
      endif
c
c      Next section computes N-body initial conditions
c
      idum = iseed
      rc3 = 2.d0*g*mtot*eps/(alpha1*(wc1 + s0*w0c1))
      rc = rc3**(1.d0/3.d0)
      rc2 = rc**2
      j22 = 2.d0*eps/(rc2*alpha1)
      pe = 0.d0
      ke = 0.d0
      sumx = 0.d0
      sumy = 0.d0
      sumz = 0.d0
      write (6,*) 'Generating N-body data...'
      do 40 i = 1,n
         if (n.le.1) goto 40
   15    continue
         rstar = rlag*ran2(idum)
         costh = 2.d0*ran2(idum)-1.d0
         phi = 2.d0*pi*ran2(idum)
         sinth = sqrt(1.d0-costh**2)
         xstar = rstar*sinth*cos(phi)
         ystar = rstar*sinth*sin(phi)
         zstar = rstar*costh
         r2s2 = (alpha3/alpha1-0.5d0)*(zstar**2-0.5d0*(xstar**2+
     &        ystar**2))+ 0.75d0*(xstar**2-ystar**2)
         r2 = xstar**2 + ystar**2 + zstar**2
         s2 = r2s2/r2
         r = sqrt(r2)
         if (r.lt.rking) then
            if (r.lt.sqrt(yking(1,1))) then
               w1 = x(1) + y(1,1)*s0 + y(3,1)*s2 - eps*r**2*
     &              (2.d0/3.d0)*(s0/3.d0 + s2)
               w = w0 - (r2/yking(1,1))*(w0 - w1)
            else
               j = 1
   20          continue
               j = j + 1
               if (r.gt.sqrt(yking(1,j))) goto 20
               wj1 = x(j-1) + y(1,j-1)*s0 + y(3,j-1)*s2 - eps*r**2*
     &              (2.d0/3.d0)*(s0/3.d0 + s2)
               wj = x(j) + y(1,j)*s0 + y(3,j)*s2 - eps*r**2*
     &              (2.d0/3.d0)*(s0/3.d0 + s2)
               w = wj1 + (r2-yking(1,j-1))*(wj-wj1)/
     &              (yking(1,j)-yking(1,j-1))
            endif
         else
            w = wc1/r + wc2
     &           +s0*(w0c1/r + w0c2) + s2*w2c1/r**3 - eps*r**2*
     &           (2.d0/3.d0)*(s0/3.d0 + s2)
         endif
         if (w.gt.0.d0) then
            rhost = densty(w)/den
         else
            rhost = 0
         endif
         if (rstar**2*rhost.lt.fmax*ran2(idum)) go to 15
c        now choose speed
         vmax = sqrt(2.d0*w)
   30    continue
         speed = vmax*ran2(idum)
         fstar = speed**2*(exp(-0.5d0*speed**2)-exp(-w))
         if (fstar.lt.2.d0*ran2(idum)/exp(1.d0)) go to 30
         costh = 2.d0*ran2(idum)-1.d0
         phi = 2.d0*pi*ran2(idum)
         sinth = sqrt(1.d0-costh**2)
         speed = speed/sqrt(j22)
         xstar = xstar*rc
         ystar = ystar*rc
         zstar = zstar*rc
         ustar = speed*sinth*cos(phi)
         vstar = speed*sinth*sin(phi)
         wstar = speed*costh
         mstar = mtot/n
         write (9,102) mstar,xstar,ystar,zstar,ustar,vstar,wstar
  102    format (1p7e10.2)
         coord(i,1) = xstar
         coord(i,2) = ystar
         coord(i,3) = zstar
         if (i.gt.1) then
            do j = 1,i-1
               r2 = (coord(j,1)-xstar)**2 + (coord(j,2)-ystar)**2 + 
     &              (coord(j,3)-zstar)**2
               pe = pe - 1.d0/sqrt(r2)
            enddo
         endif
         ke = ke + mstar*speed**2
         sumx = sumx + mstar*xstar**2
         sumy = sumy + mstar*ystar**2
         sumz = sumz + mstar*zstar**2
   40 continue
      pe = g*pe*mstar**2
      ke = ke*0.5d0
      xmax = rc*rlag
      ymax = rc*ymax
      zmax = rc*zmax
      write (6,*) 'radii (user''s units):',xmax,ymax,zmax
      if (n.gt.1) then
         rvir = -g*mtot**2/(2.d0*pe)
         write (6,*) 'Virial radius (user''s units):',rvir
c
c    Following section assumes units of solar masses, km/sec and pc,
c    and A = 14.4, B = -12.0 km/sec/kpc
c
         omega = -0.0264
         kevir = -0.5d0*pe+ 0.5d0*alpha1*sumx + 0.5d0*alpha3*sumz +
     &     omega**2*(sumx + sumy)
         write (6,*) 'Kinetic energy:',ke,' potential energy',pe
         write (6,*) 'Kinetic energy from virial theorem (assuming ',
     &     'certain units):',kevir
      endif

      stop
      END
      SUBROUTINE DIFEQ(K,K1,K2,JSF,IS1,ISF,INDEXV,NE,S,NSI,NSJ,Y,NYJ,NYK
     *)
	implicit double precision (a-h,o-z)
      PARAMETER(M=51)
      COMMON X(M),H,MM,N,C2,ANORM
	common /params/ w0,alpha1,alpha3,rho(m),rhodot(m),rlag,
     &     yking(2,m),den,rking
      DIMENSION Y(NYJ,NYK),S(NSI,NSJ),INDEXV(NYJ)
      g2(i) = -2.25d0*yking(2,i)**2*(rhodot(i)*(y(1,i) + (2.d0/9.d0)*
     &     exp(y(5,i))*yking(1,i)) - y(2,i)*rho(i))/yking(1,i)
      g4(i) = 0.25d0*yking(2,i)**2*(6.d0*y(3,i)/yking(1,i)-
     &     9.d0*rhodot(i)*(y(3,i) + (2.d0/3.d0)*exp(y(5,i))*yking(1,i))
     &		+ 9.d0*y(4,i)*rho(i))/yking(1,i)
      IF(K.EQ.K1) THEN
	do i = 3,5
		do j = 6,10
			s(i,j) = 0.d0
		enddo
	enddo
	eps = -exp(y(5,1))
	b0 = 0.1d0*eps*rhodot(1)
	s(3,5+1) = 1.d0
	s(3,5+5) = -b0*yking(1,1)**2
	s(3,jsf) = y(1,1) - b0*yking(1,1)**2
	s(4,5+2) = 1.d0
	s(4,5+5) = -2.d0*b0*yking(1,1)*yking(2,1)
	s(4,jsf) = y(2,1) - 2.d0*b0*yking(1,1)*yking(2,1)
	s(5,5+3) = yking(2,1)
	s(5,5+4) = -yking(1,1)
	s(5,jsf) = y(3,1)*yking(2,1) - yking(1,1)*y(4,1)
      ELSE IF(K.GT.K2) THEN
	do i = 1,2
		do j = 6,10
			s(i,j) = 0.d0
		enddo
	enddo
	s(1,5+3) = 1.5d0*yking(2,m)
	s(1,5+4) = yking(1,m)
	s(1,jsf) = yking(1,m)*y(4,m) + 1.5d0*y(3,m)*yking(2,m)
c	now find radius where acceleration vanishes
	r = sqrt(yking(1,m))
	s0 = 1.5d0*(1.d0+alpha3/alpha1)
	s2 = -0.5d0*alpha3/alpha1 + 1.d0
	iter = 0
        eps = -exp(y(5,m))
	wc1 = -2.d0*yking(1,m)**1.5d0/yking(2,m)
	wc2 = -wc1/sqrt(yking(1,m))
	w0c1 = wc1*y(2,m)
	w0c2 = y(1,m) - w0c1/sqrt(yking(1,m))
	w2c1 = yking(1,m)**1.5d0*y(3,m)
10	continue
	iter = iter + 1
	if (iter.gt.20) then
           write (6,*) (y(i,m),i=1,5)
           write (6,*) 'too many iterations finding r'
           stop
        endif
	f = -wc1/r**2 +s0*(-w0c1/r**2)
     &		-s2*3.d0*w2c1/r**4 - 2.d0*eps*r
	fdash = 2.d0*wc1/r**3 + 2.d0*s0*w0c1
     &		/r**3 +12.d0*s2*w2c1/r**5 -
     &		2.d0*eps
	rnew = r - f/fdash
	if (abs((r-rnew)/rnew).gt.1.d-4) then
		r = rnew
		go to 10
	endif
	r = rnew
	rlag = r
        temp = 1.d0 - sqrt(yking(1,m))/r
	s(2,5+1) = s0
	s(2,5+2) = s0*2.d0*(yking(1,m)/yking(2,m))*temp
	s(2,5+3) = s2*yking(1,m)**1.5d0/r**3
	s(2,5+5) = exp(y(5,m))*r**2
	s(2,jsf) = 2.d0*yking(1,m)*temp/yking(2,m) +
     &		s0*(y(1,m) + 2.d0*yking(1,m)*y(2,m)*temp/yking(2,m))
     &		+ s2*yking(1,m)**1.5d0*y(3,m)/r**3
     &		- eps*r**2
      ELSE
	do i = 1,5
		do j = 1,10
			s(i,j)=0.d0
		enddo
	enddo
	halfh = 0.5d0*(x(k)-x(k-1))
	s(1,1) = -1.d0
	s(1,2) = -halfh
	s(1,5+1)=1.d0
	s(1,5+2)=-halfh
	s(1,jsf)=y(1,k)-y(1,k-1)-halfh*(y(2,k)+y(2,k-1))
  	s(2,1) = -halfh*(-2.25d0*(yking(2,k-1)**2/yking(1,k-1))*
     &       rhodot(k-1))
	s(2,2) = -1.d0 - halfh*(-2.25d0*(yking(2,k-1)**2/yking(1,k-1))*
     &       (-rho(k-1)))
	s(2,5) = -halfh*(-2.25d0*(yking(2,k-1)**2/yking(1,k-1))*
     &       rhodot(k-1)*
     &       (2.d0/9.d0)*exp(y(5,k-1))*yking(1,k-1))
  	s(2,5+1) = -halfh*(-2.25d0*(yking(2,k)**2/yking(1,k))*rhodot(k))
	s(2,5+2) = 1.d0 - halfh*(-2.25d0*(yking(2,k)**2/yking(1,k))*
     &       (-rho(k)))
	s(2,5+5) = -halfh*(-2.25d0*(yking(2,k)**2/yking(1,k))*rhodot(k)*
     &       (2.d0/9.d0)*exp(y(5,k))*yking(1,k))
	s(2,jsf) = y(2,k)-y(2,k-1)-halfh*(g2(k-1)+g2(k))
	s(3,3) = -1.d0
	s(3,4) = -halfh
	s(3,5+3)=1.d0
	s(3,5+4)=-halfh
	s(3,jsf)=y(3,k)-y(3,k-1)-halfh*(y(4,k)+y(4,k-1))
  	s(4,3) = -halfh*(0.25d0*(yking(2,k-1)**2/yking(1,k-1))*(6.d0/
     &       yking(1,k-1) 
     &       -9.d0*rhodot(k-1)))
	s(4,4) = -1.d0 - halfh*(0.25d0*(yking(2,k-1)**2/yking(1,k-1))*
     &       9.d0*
     &       rho(k-1))
	s(4,5) = -halfh*(0.25d0*(yking(2,k-1)**2/yking(1,k-1))*(-9.d0)*
     &       rhodot(k-1)*(2.d0/3.d0)*exp(y(5,k-1))*yking(1,k-1))
  	s(4,5+3) = -halfh*(0.25d0*(yking(2,k)**2/yking(1,k))*(6.d0/
     &       yking(1,k) 
     &       -9.d0*rhodot(k)))
	s(4,5+4) = 1.d0 - halfh*(0.25d0*(yking(2,k)**2/yking(1,k))*9.d0*
     &       rho(k))
	s(4,5+5) = -halfh*(0.25d0*(yking(2,k)**2/yking(1,k))*(-9.d0)*
     &       rhodot(k)*(2.d0/3.d0)*exp(y(5,k))*yking(1,k))
	s(4,jsf) = y(4,k)-y(4,k-1)-halfh*(g4(k-1)+g4(k))
	s(5,5) = -1.d0
	s(5,5+5) = 1.d0
	s(5,jsf) = -y(5,k-1) + y(5,k)
      ENDIF
      RETURN
      END
      SUBROUTINE SOLVDE(ITMAX,CONV,SLOWC,SCALV,INDEXV,NE,NB,M,
     *    Y,NYJ,NYK,C,NCI,NCJ,NCK,S,NSI,NSJ)
	implicit double precision (a-h,o-z)
      PARAMETER (NMAX=10)
      DIMENSION Y(NYJ,NYK),C(NCI,NCJ,NCK),S(NSI,NSJ),SCALV(NYJ),INDEXV(N
     *YJ)
      DIMENSION ERMAX(NMAX),KMAX(NMAX)
      K1=1
      K2=M
      NVARS=NE*M
      J1=1
      J2=NB
      J3=NB+1
      J4=NE
      J5=J4+J1
      J6=J4+J2
      J7=J4+J3
      J8=J4+J4
      J9=J8+J1
      IC1=1
      IC2=NE-NB
      IC3=IC2+1
      IC4=NE
      JC1=1
      JCF=IC3
      DO 16 IT=1,ITMAX
        K=K1
        CALL DIFEQ(K,K1,K2,J9,IC3,IC4,INDEXV,NE,S,NSI,NSJ,Y,NYJ,NYK)
        CALL PINVS(IC3,IC4,J5,J9,JC1,K1,C,NCI,NCJ,NCK,S,NSI,NSJ)
        DO 11 K=K1+1,K2
          KP=K-1
          CALL DIFEQ(K,K1,K2,J9,IC1,IC4,INDEXV,NE,S,NSI,NSJ,Y,NYJ,NYK)
          CALL RED(IC1,IC4,J1,J2,J3,J4,J9,IC3,JC1,JCF,KP,
     *        C,NCI,NCJ,NCK,S,NSI,NSJ)
          CALL PINVS(IC1,IC4,J3,J9,JC1,K,C,NCI,NCJ,NCK,S,NSI,NSJ)
11      CONTINUE
        K=K2+1
        CALL DIFEQ(K,K1,K2,J9,IC1,IC2,INDEXV,NE,S,NSI,NSJ,Y,NYJ,NYK)
        CALL RED(IC1,IC2,J5,J6,J7,J8,J9,IC3,JC1,JCF,K2,
     *      C,NCI,NCJ,NCK,S,NSI,NSJ)
        CALL PINVS(IC1,IC2,J7,J9,JCF,K2+1,C,NCI,NCJ,NCK,S,NSI,NSJ)
        CALL BKSUB(NE,NB,JCF,K1,K2,C,NCI,NCJ,NCK)
        ERR=0.
        DO 13 J=1,NE
          JV=INDEXV(J)
          ERMAX(J)=0.
          ERRJ=0.
          KMAX(J)=0
          VMAX=0.
          DO 12 K=K1,K2
            VZ=ABS(C(J,1,K))
            IF(VZ.GT.VMAX) THEN
               VMAX=VZ
               KM=K
            ENDIF
            ERRJ=ERRJ+VZ
12        CONTINUE
          ERR=ERR+ERRJ/SCALV(JV)
          ERMAX(J)=C(J,1,KM)/SCALV(JV)
          KMAX(J)=KM
13      CONTINUE
        ERR=ERR/NVARS
        FAC=SLOWC/MAX(SLOWC,ERR)
        DO 15 JV=1,NE
          J=INDEXV(JV)
          DO 14 K=K1,K2
            Y(J,K)=Y(J,K)-FAC*C(JV,1,K)
14        CONTINUE
15      CONTINUE
        WRITE(6,100) IT,ERR,FAC,(KMAX(J),ERMAX(J),J=1,NE)
        IF(ERR.LT.CONV) RETURN
16    CONTINUE
      write (6,*) 'ITMAX exceeded'
      stop
100   FORMAT(1X,I4,2F12.6,(/5X,I5,F12.6))
      END
      SUBROUTINE BKSUB(NE,NB,JF,K1,K2,C,NCI,NCJ,NCK)
	implicit double precision (a-h,o-z)
      DIMENSION C(NCI,NCJ,NCK)
      NBF=NE-NB
      DO 13 K=K2,K1,-1
        KP=K+1
        DO 12 J=1,NBF
          XX=C(J,JF,KP)
          DO 11 I=1,NE
            C(I,JF,K)=C(I,JF,K)-C(I,J,K)*XX
11        CONTINUE
12      CONTINUE
13    CONTINUE
      DO 16 K=K1,K2
        KP=K+1
        DO 14 I=1,NB
          C(I,1,K)=C(I+NBF,JF,K)
14      CONTINUE
        DO 15 I=1,NBF
          C(I+NB,1,K)=C(I,JF,KP)
15      CONTINUE
16    CONTINUE
      RETURN
      END
      SUBROUTINE PINVS(IE1,IE2,JE1,JSF,JC1,K,C,NCI,NCJ,NCK,S,NSI,NSJ)
	implicit double precision (a-h,o-z)
      PARAMETER (ZERO=0.,ONE=1.,NMAX=10)
      DIMENSION C(NCI,NCJ,NCK),S(NSI,NSJ),PSCL(NMAX),INDXR(NMAX)
      JE2=JE1+IE2-IE1
      JS1=JE2+1
      DO 12 I=IE1,IE2
        BIG=ZERO
        DO 11 J=JE1,JE2
          IF(ABS(S(I,J)).GT.BIG) BIG=ABS(S(I,J))
11      CONTINUE
        IF(BIG.EQ.ZERO) then
           write(6,*) 'Singular matrix, row all 0'
           stop
        endif
        PSCL(I)=ONE/BIG
        INDXR(I)=0
12    CONTINUE
      DO 18 ID=IE1,IE2
        PIV=ZERO
        DO 14 I=IE1,IE2
          IF(INDXR(I).EQ.0) THEN
            BIG=ZERO
            DO 13 J=JE1,JE2
              IF(ABS(S(I,J)).GT.BIG) THEN
                JP=J
                BIG=ABS(S(I,J))
              ENDIF
13          CONTINUE
            IF(BIG*PSCL(I).GT.PIV) THEN
              IPIV=I
              JPIV=JP
              PIV=BIG*PSCL(I)
            ENDIF
          ENDIF
14      CONTINUE
        IF(S(IPIV,JPIV).EQ.ZERO) then
           write (6,*) 'Singular matrix'
           stop
        endif
        INDXR(IPIV)=JPIV
        PIVINV=ONE/S(IPIV,JPIV)
        DO 15 J=JE1,JSF
          S(IPIV,J)=S(IPIV,J)*PIVINV
15      CONTINUE
        S(IPIV,JPIV)=ONE
        DO 17 I=IE1,IE2
          IF(INDXR(I).NE.JPIV) THEN
            IF(S(I,JPIV).NE.ZERO) THEN
              DUM=S(I,JPIV)
              DO 16 J=JE1,JSF
                S(I,J)=S(I,J)-DUM*S(IPIV,J)
16            CONTINUE
              S(I,JPIV)=ZERO
            ENDIF
          ENDIF
17      CONTINUE
18    CONTINUE
      JCOFF=JC1-JS1
      ICOFF=IE1-JE1
      DO 21 I=IE1,IE2
        IROW=INDXR(I)+ICOFF
        DO 19 J=JS1,JSF
          C(IROW,J+JCOFF,K)=S(I,J)
19      CONTINUE
21    CONTINUE
      RETURN
      END
      SUBROUTINE RED(IZ1,IZ2,JZ1,JZ2,JM1,JM2,JMF,IC1,JC1,JCF,KC,
     *    C,NCI,NCJ,NCK,S,NSI,NSJ)
	implicit double precision (a-h,o-z)
      DIMENSION C(NCI,NCJ,NCK),S(NSI,NSJ)
      LOFF=JC1-JM1
      IC=IC1
      DO 14 J=JZ1,JZ2
        DO 12 L=JM1,JM2
          VX=C(IC,L+LOFF,KC)
          DO 11 I=IZ1,IZ2
            S(I,L)=S(I,L)-S(I,J)*VX
11        CONTINUE
12      CONTINUE
        VX=C(IC,JCF,KC)
        DO 13 I=IZ1,IZ2
          S(I,JMF)=S(I,JMF)-S(I,J)*VX
13      CONTINUE
        IC=IC+1
14    CONTINUE
      RETURN
      END
c
c
      double precision FUNCTION ERF(X)
      implicit double precision (a-h,o-z)
	half = 0.5d0
      IF(X.LT.0.)THEN
        ERF=-GAMMP(half,X**2)
      ELSE
        ERF=GAMMP(half,X**2)
      ENDIF
      RETURN
      END
c
c
      double precision FUNCTION GAMMP(A,X)
      implicit double precision (a-h,o-z)
      IF(X.LT.0..OR.A.LE.0.) stop
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMMP,A,X,GLN)
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMP=1.-GAMMCF
      ENDIF
      RETURN
      END
c
c
      SUBROUTINE GSER(GAMSER,A,X,GLN)
      implicit double precision (a-h,o-z)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
        IF(X.LT.0.) stop
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    CONTINUE
      write (6,*) 'A too large, ITMAX too small'
      stop
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END
c
c
      double precision FUNCTION GAMMLN(XX)
      implicit double precision (a-h,o-z)
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
c
c
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      implicit double precision (a-h,o-z)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO 11 N=1,ITMAX
        AN=FLOAT(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.0.)THEN
          FAC=1./A1
          G=B1*FAC
          IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      write (6,*) 'A too large, ITMAX too small'
      stop
1     GAMMCF=EXP(-X+A*LOG(X)-GLN)*G
      RETURN
      END
c
      FUNCTION RAN2(IDUM)
      PARAMETER (M=714025,IA=1366,IC=150889,RM=1.4005112E-6)
      implicit double precision (a-h,o-z)
      DIMENSION IR(97)
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        IDUM=MOD(IC-IDUM,M)
        DO 11 J=1,97
          IDUM=MOD(IA*IDUM+IC,M)
          IR(J)=IDUM
11      CONTINUE
        IDUM=MOD(IA*IDUM+IC,M)
        IY=IDUM
      ENDIF
      J=1+(97*IY)/M
      IF(J.GT.97.OR.J.LT.1) stop
      IY=IR(J)
      RAN2=IY*RM
      IDUM=MOD(IA*IDUM+IC,M)
      IR(J)=IDUM
      RETURN
      END
c
c
      SUBROUTINE ODEINT(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,DERIVS,RK
     *QC)
      PARAMETER (MAXSTP=10000,NMAX=10,TWO=2.0,ZERO=0.0,TINY=1.E-30)
	implicit double precision(a-h,o-z)
      COMMON /PATH/ KMAX,KOUNT,DXSAV,XP(200),YP(10,200)
      DIMENSION YSTART(NVAR),YSCAL(NMAX),Y(NMAX),DYDX(NMAX)
      X=X1
      H=SIGN(H1,X2-X1)
      NOK=0
      NBAD=0
      KOUNT=0
      DO 11 I=1,NVAR
        Y(I)=YSTART(I)
11    CONTINUE
      XSAV=X-DXSAV*TWO
      DO 16 NSTP=1,MAXSTP
        CALL DERIVS(X,Y,DYDX)
        DO 12 I=1,NVAR
          YSCAL(I)=ABS(Y(I))+ABS(H*DYDX(I))+TINY
12      CONTINUE
        IF(KMAX.GT.0)THEN
          IF(ABS(X-XSAV).GT.ABS(DXSAV)) THEN
            IF(KOUNT.LT.KMAX-1)THEN
              KOUNT=KOUNT+1
              XP(KOUNT)=X
              DO 13 I=1,NVAR
                YP(I,KOUNT)=Y(I)
13            CONTINUE
              XSAV=X
            ENDIF
          ENDIF
        ENDIF
        IF((X+H-X2)*(X+H-X1).GT.ZERO) H=X2-X
        CALL RKQC(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
        IF(HDID.EQ.H)THEN
          NOK=NOK+1
        ELSE
          NBAD=NBAD+1
        ENDIF
        IF((X-X2)*(X2-X1).GE.ZERO)THEN
          DO 14 I=1,NVAR
            YSTART(I)=Y(I)
14        CONTINUE
          IF(KMAX.NE.0)THEN
            KOUNT=KOUNT+1
            XP(KOUNT)=X
            DO 15 I=1,NVAR
              YP(I,KOUNT)=Y(I)
15          CONTINUE
          ENDIF
          RETURN
        ENDIF
        IF(ABS(HNEXT).LT.HMIN) then
           write (6,*) 'Stepsize smaller than minimum.'
           stop
        endif
        H=HNEXT
16    CONTINUE
      write (6,*) 'Too many steps.'
      stop
      END
c
c
      SUBROUTINE RKQC(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
      PARAMETER (NMAX=10,FCOR=.0666666667,
     *    ONE=1.,SAFETY=0.9,ERRCON=6.E-4)
	implicit double precision (a-h,o-z)
      EXTERNAL DERIVS
      DIMENSION Y(N),DYDX(N),YSCAL(N),YTEMP(NMAX),YSAV(NMAX),DYSAV(NMAX)
      PGROW=-0.20
      PSHRNK=-0.25
      XSAV=X
      DO 11 I=1,N
        YSAV(I)=Y(I)
        DYSAV(I)=DYDX(I)
11    CONTINUE
      H=HTRY
1     HH=0.5*H
      CALL RK4(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)
      X=XSAV+HH
      CALL DERIVS(X,YTEMP,DYDX)
      CALL RK4(YTEMP,DYDX,N,X,HH,Y,DERIVS)
      X=XSAV+H
      IF(X.EQ.XSAV) then
         write (6,*) 'Stepsize not significant in RKQC.'
         stop
      endif
      CALL RK4(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
      ERRMAX=0.
      DO 12 I=1,N
        YTEMP(I)=Y(I)-YTEMP(I)
        ERRMAX=MAX(ERRMAX,ABS(YTEMP(I)/YSCAL(I)))
12    CONTINUE
      ERRMAX=ERRMAX/EPS
      IF(ERRMAX.GT.ONE) THEN
        H=SAFETY*H*(ERRMAX**PSHRNK)
        GOTO 1
      ELSE
        HDID=H
        IF(ERRMAX.GT.ERRCON)THEN
          HNEXT=SAFETY*H*(ERRMAX**PGROW)
        ELSE
          HNEXT=4.*H
        ENDIF
      ENDIF
      DO 13 I=1,N
        Y(I)=Y(I)+YTEMP(I)*FCOR
13    CONTINUE
      RETURN
      END
c
c
      SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT,DERIVS)
      PARAMETER (NMAX=10)
	implicit double precision (a-h,o-z)
      DIMENSION Y(N),DYDX(N),YOUT(N),YT(NMAX),DYT(NMAX),DYM(NMAX)
      HH=H*0.5
      H6=H/6.
      XH=X+HH
      DO 11 I=1,N
        YT(I)=Y(I)+HH*DYDX(I)
11    CONTINUE
      CALL DERIVS(XH,YT,DYT)
      DO 12 I=1,N
        YT(I)=Y(I)+HH*DYT(I)
12    CONTINUE
      CALL DERIVS(XH,YT,DYM)
      DO 13 I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
13    CONTINUE
      CALL DERIVS(X+H,YT,DYT)
      DO 14 I=1,N
        YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
14    CONTINUE
      RETURN
      END
c
c
	subroutine derivs (x,y,dydx)
	implicit double precision (a-h,o-z)
	parameter (nvar=2,m=51)
      common /params/ w0,alpha1,alpha3,rho(m),rhodot(m),rlag,yking(2,m),
     &     den,rking
	dimension y(nvar),dydx(nvar)
	pi = 4.d0*datan(1.d0)
	if (x.ge.0.d0) then
	      rhox =-sqrt(x)*(x+1.5d0)+0.75*sqrt(pi)*exp(x)*erf(sqrt(x))
	else
		rhox = 0.d0
	endif
	dydx(1) = y(2)
	dydx(2) = 0.25d0*y(2)**2*(6.d0+9.d0*y(2)*rhox/den)/y(1)	
	return
	end
c
c
      double precision function sigma(xstar,zstar,y)
        implicit double precision (a-h,o-z)
        parameter (ne=5,m=51)
      common /params/ w0,alpha1,alpha3,rho(m),rhodot(m),rlag,yking(2,m),
     &     den,rking
      COMMON X(M),H,MM,N,C2,ANORM
      dimension y(ne,m)
      densty(z) = -sqrt(z)*(z+1.5d0)+0.75*sqrt(pi)*exp(z)*erf(sqrt(z))
        eps = -exp(y(5,m))
      s0 = 1.5d0*(1.d0+alpha3/alpha1)
	wc1 = -2.d0*yking(1,m)**1.5d0/yking(2,m)
	wc2 = -wc1/sqrt(yking(1,m))
	w0c1 = wc1*y(2,m)
	w0c2 = y(1,m) - w0c1/sqrt(yking(1,m))
	w2c1 = yking(1,m)**1.5d0*y(3,m)
        pi = 4.d0*atan(1.d0)
        sum = 0.d0
        ystar = 0.d0
        coeff = 0.5d0
   10   continue
         r2s2 = (alpha3/alpha1-0.5d0)*(zstar**2-0.5d0*(xstar**2+
     &        ystar**2))+ 0.75d0*(xstar**2-ystar**2)
         r2 = xstar**2 + ystar**2 + zstar**2
         if (r2.gt.0.d0) then
            s2 = r2s2/r2
         else
            s2 = 1.d0
         endif
         r = sqrt(r2)
         if (r.lt.rking) then
            if (r.lt.sqrt(yking(1,1))) then
               w1 = x(1) + y(1,1)*s0 + y(3,1)*s2 - eps*r**2*
     &              (2.d0/3.d0)*(s0/3.d0 + s2)
               w = w0 - (r2/yking(1,1))*(w0 - w1)
            else
               j = 1
   20          continue
               j = j + 1
               if (r.gt.sqrt(yking(1,j))) goto 20
               wj1 = x(j-1) + y(1,j-1)*s0 + y(3,j-1)*s2 - eps*r**2*
     &              (2.d0/3.d0)*(s0/3.d0 + s2)
               wj = x(j) + y(1,j)*s0 + y(3,j)*s2 - eps*r**2*
     &              (2.d0/3.d0)*(s0/3.d0 + s2)
               w = wj1 + (r2-yking(1,j-1))*(wj-wj1)/
     &              (yking(1,j)-yking(1,j-1))
            endif
         else
            w = wc1/r + wc2
     &           +s0*(w0c1/r + w0c2) + s2*w2c1/r**3 - eps*r**2*
     &           (2.d0/3.d0)*(s0/3.d0 + s2)
         endif
         if (w.gt.0.d0) then
            rhost = densty(w)/den
            sum = sum + coeff*rhost
            if (coeff.lt.1.d0) coeff=1.d0
            ystar = ystar + 0.1
            goto 10
         else
            sigma = sum
            return
         endif
         end
