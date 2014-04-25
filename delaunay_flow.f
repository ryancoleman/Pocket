c	Delaunay_flow.f		Version 1 11/24/2000	Patrice Koehl
c
c	This subroutine defines the flow of each tetrahedron of the
c	Delaunay
c
	subroutine delaunay_flow
c
	include 'pocket.h'
c
	integer i,j,k,l,m,n,p
	integer	ntetra
	integer idx
	integer	iflow(4)
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
c	Information on the tetrahedra of the regular
c	triangulation
c
	integer tetra(4,ntetra_max)
	integer tetra_neighbour(4,ntetra_max)
c
	integer*1 tetra_flow(4,ntetra_max)
c
	real*8	ra,rb,rc,rd,wa,wb,wc,wd
	real*8	scale,eps
	real*8	a(3),b(3),c(3),d(3)
	real*8	radius(npointmax),weight(npointmax)
	real*8	coord(3*npointmax)
c
	common  /xyz_vertex/	coord,radius,weight
	common  /tetra_zone/    ntetra,tetra,tetra_neighbour
	common  /tetra_stat/    tetra_status,tetra_orient,
     1                          tetra_nindex
	common  /tetra_flux/	tetra_flow
	common  /gmp_info/	scale,eps
c
c	Loop over all tetrahedra:
c
	do 400 idx = 1,ntetra
c
c		"Dead" tetrahedron are ignored
c
		if(tetra_status(idx).eq.0) goto 400
c
		i = tetra(1,idx)
		j = tetra(2,idx)
		k = tetra(3,idx)
		l = tetra(4,idx)
c
		do 100 m = 1,3
			a(m) = coord(3*(i-1)+m)
			b(m) = coord(3*(j-1)+m)
			c(m) = coord(3*(k-1)+m)
			d(m) = coord(3*(l-1)+m)
100		continue
c
		ra = radius(i)
		rb = radius(j)
		rc = radius(k)
		rd = radius(l)
c
		wa = weight(i)
		wb = weight(j)
		wc = weight(k)
		wd = weight(l)
c
		do 200 m = 1,4
			iflow(m) = 0
			n = tetra_neighbour(m,idx)
			p = tetra_nindex(m,idx)
			if(n.ne.0.and.n.lt.idx) then
				if(tetra_flow(p,n).eq.1) then
					iflow(m) = -1
				endif
			endif
200		continue
c
		call flow(a,b,c,d,ra,rb,rc,rd,wa,wb,wc,wd,
     1		iflow,eps,scale)
c
		do 300 m = 1,4
			if(iflow(m).eq.1) then
				tetra_flow(m,idx) = 1
			else
				tetra_flow(m,idx) = 0
			endif
300		continue
c
400	continue
c
	return
	end
