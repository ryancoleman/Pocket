c	Constants.def		Version 1 8/29/2000		Patrice Koehl
c
c	This file contains the definition of some constants
c
	real*8	pi,twopi,fourpi,precision
	integer	ibuild
c
	save	ibuild,pi,twopi,fourpi,precision
c
	data ibuild /0/
c
	if(ibuild.eq.0) then
c
		pi = acos(-1.d0)
		twopi = 2*pi
		fourpi = 4*pi
		precision = 1.d-8
		ibuild = 1
c
	endif
c
