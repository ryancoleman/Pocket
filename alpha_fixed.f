c	Alpha_fixed.f		Version 1 11/24/2000	Patrice Koehl
c
c	This file contains a suite of routine for generating the alpha
c	shape based on a given weighted Delaunay triangulation, and
c	a fixed value of Alpha (0 when we are interested in a union of balls)
c
c	Copyright (C) 2002 Patrice Koehl
c
c	This library is free software; you can redistribute it and/or
c	modify it under the terms of the GNU Lesser General Public
c	License as published by the Free Software Foundation; either
c	version 2.1 of the License, or (at your option) any later version.
c
c	This library is distributed in the hope that it will be useful,
c	but WITHOUT ANY WARRANTY; without even the implied warranty of
c	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
c	Lesser General Public License for more details.
c
c	You should have received a copy of the GNU Lesser General Public
c	License along with this library; if not, write to the Free Software
c	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
c
c
c	Alpha_fixed.f		Version 1 12/1/2000	Patrice Koehl
c
c	This subroutine builds the alpha complex based on the weighted
c	Delaunay triangulation
c
	subroutine alpha_fixed(alpha)
c
	include 'pocket.h'
c
	integer i,j,k,l,m,f,e
	integer	ntetra,ntrig,nedge
	integer	npoints,nvertex
	integer idx
	integer	pair(6),triangle(4)
	integer	hullinfo(4)
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
c	Information on the tetrahedra of the regular
c	triangulation
c
	integer tetra(4,ntetra_max)
	integer tetra_neighbour(4,ntetra_max)
	integer tetra_link(4,ntetra_max)
c
c	Information on the triangles of the regular 
c	triangulation
c
	integer trig(3,ntrig_max)
	integer trig_link(3,ntrig_max)
	integer*1 trig_status(ntrig_max)
	integer*1 trig_coef(ntrig_max)
c
c	Information on the edges of the regular
c	triangulation
c
	integer edge(2,nedge_max)
	integer*1 edge_status(nedge_max)
c
c	Information on the vertices
c
	integer redinfo(npointmax)
c
c	Information on convex hull
c
	integer*1 tetra_hull(ntetra_max),trig_hull(ntrig_max)
	integer*1 edge_hull(nedge_max),vertex_hull(npointmax)
c
	real*8	ra,rb,rc,rd,wa,wb,wc,wd
	real*8	alpha,scale,eps
	real*8	a(3),b(3),c(3),d(3)
	real*8	radius(npointmax),weight(npointmax)
	real*8	coord(3*npointmax)
c
	logical r_tetra(ntetra_max),r_face(ntrig_max)
	logical a_tetra(ntetra_max),a_face(ntrig_max)
	logical r_edge(nedge_max),r_vertex(npointmax)
	logical a_edge(nedge_max),a_vertex(npointmax)
	logical testr(15),testa(15)
c
	common  /xyz_vertex/	coord,radius,weight
	common  /tetra_zone/    ntetra,tetra,tetra_neighbour
	common  /tetra_stat/    tetra_status,tetra_orient,
     1                          tetra_nindex
	common  /trig_zone/     ntrig,trig
	common  /trig_stat/	trig_status,trig_coef
	common  /edge_zone/     nedge,edge
	common  /edge_stat/	edge_status
	common  /vertex_zone/   npoints,nvertex,redinfo
	common  /links/         tetra_link,trig_link
	common  /gmp_info/	scale,eps
	common  /hull/		tetra_hull,trig_hull,edge_hull,
     1				vertex_hull
c
c	Initialise the two arrays r_tetra and a_tetra:
c
c	r_tetra(i) is .true. if the simplex i is in the
c	Alpha complex, .false. otherwise
c	a_tetra(i) is .true. if the simple i is "attached"
c	.false. otherwise (see alf_tetra for details)
c
	do 200 i = 1,ntetra
		r_tetra(i) = .False.
		a_tetra(i) = .False.
200	continue
	do 300 i = 1,ntrig
		r_face(i) = .False.
		a_face(i) = .False.
300	continue
	do 400 i = 1,nedge
		r_edge(i) = .False.
		a_edge(i) = .False.
400	continue
	do 500 i = 1,npoints
		r_vertex(i) = .False.
		a_vertex(i) = .False.
500	continue
c
c	Loop over all tetrahedra:
c
	do 1100 idx = 1,ntetra
c
c		"Dead" tetrahedron are ignored
c
		if(tetra_status(idx).eq.0) goto 1100
c
		i = tetra(1,idx)
		j = tetra(2,idx)
		k = tetra(3,idx)
		l = tetra(4,idx)
c
		do 600 m = 1,3
			a(m) = coord(3*(i-1)+m)
			b(m) = coord(3*(j-1)+m)
			c(m) = coord(3*(k-1)+m)
			d(m) = coord(3*(l-1)+m)
600		continue
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
c	Get triangles (ijk), (ijl), (ikl) and (jkl)
c	(in that order)
c
		triangle(1) = tetra_link(4,idx)
		triangle(2) = tetra_link(3,idx)
		triangle(3) = tetra_link(2,idx)
		triangle(4) = tetra_link(1,idx)
c
		hullinfo(1) = trig_hull(triangle(1))
		hullinfo(2) = trig_hull(triangle(2))
		hullinfo(3) = trig_hull(triangle(3))
		hullinfo(4) = trig_hull(triangle(4))
c
c	Get edges (ij),(ik),(il),(jk),(jl),(kl)
c	(in that order)
c
		pair(1) = trig_link(3,triangle(1))
		pair(2) = trig_link(2,triangle(1))
		pair(3) = trig_link(2,triangle(2))
		pair(4) = trig_link(1,triangle(1))
		pair(5) = trig_link(1,triangle(2))
		pair(6) = trig_link(1,triangle(3))
c
		testr(1) = r_tetra(idx)
		testa(1) = a_tetra(idx)
c
		do 700 m = 1,4
			testr(1+m) = r_face(triangle(m))
			testa(1+m) = a_face(triangle(m))
700		continue
c
		do 800 m = 1,6
			testr(5+m) = r_edge(pair(m))
			testa(5+m) = a_edge(pair(m))
800		continue
c
		testr(12) = r_vertex(i)
		testa(12) = a_vertex(i)
		testr(13) = r_vertex(j)
		testa(13) = a_vertex(j)
		testr(14) = r_vertex(k)
		testa(14) = a_vertex(k)
		testr(15) = r_vertex(l)
		testa(15) = a_vertex(l)
c
		call alf_tetra(a,b,c,d,ra,rb,rc,rd,wa,wb,wc,wd,
     1		testr,testa,hullinfo,eps,scale,alpha)
c
		r_tetra(idx) = testr(1)
		a_tetra(idx) = testa(1)
c
		do 900 m = 1,4
			f = triangle(m)
			r_face(f) = testr(1+m)
			a_face(f) = testa(1+m)
900		continue
c
		do 1000 m = 1,6
			e = pair(m)
			r_edge(e) = testr(5+m)
			a_edge(e) = testa(5+m)
1000		continue
c
		r_vertex(i) = testr(12)
		a_vertex(i) = testa(12)
		r_vertex(j) = testr(13)
		a_vertex(j) = testa(13)
		r_vertex(k) = testr(14)
		a_vertex(k) = testa(14)
		r_vertex(l) = testr(15)
		a_vertex(l) = testa(15)
c
1100	continue
c
	do 1200 idx = 1,npoints
		if(r_vertex(idx)) then
			redinfo(idx) = 0
		else
			redinfo(idx) = 1
		endif
1200	continue
c
	do 1300 idx = 1,nedge
		if(r_edge(idx)) then
			edge_status(idx) = 1
		else
			edge_status(idx) = 0
		endif
1300	continue
c
	do 1400 idx = 1,ntrig
		if(r_face(idx)) then
			trig_status(idx) = 1
			trig_coef(idx) = 2
		else
			trig_status(idx) = 0
			trig_coef(idx) = 0
		endif
1400	continue
c
	do 1600 idx = 1,ntetra
		if(r_tetra(idx)) then
			triangle(1) = tetra_link(1,idx)
			triangle(2) = tetra_link(2,idx)
			triangle(3) = tetra_link(3,idx)
			triangle(4) = tetra_link(4,idx)
			do 1500 m = 1,4
				trig_coef(triangle(m)) = 
     1					trig_coef(triangle(m)) -1
1500			continue
		else
			tetra_status(idx) = -tetra_status(idx)
		endif
1600	continue
c
	return
	end
c
c	alf_tetra.f	Version 1 11/24/2000	Patrice Koehl
c
c	This subroutine computes the radius R of the circumsphere containing
c	a tetrahedron [A,B,C,D], as well as check if any fourth point L 
c	(in A, B, C, D) of the tetrahedron is "hidden" by its opposite 
c	face [I,J,K], i.e. is interior to the cicumsphere of [I,J,K]
c
c	Since we are only interested at how R compares to Alpha, we do not
c	output R, rather the result of the comparison
c
c	Computation extends to all four faces of the tetrehedron, as
c	well as to all six edges of the tetrahedron
c
c	Computation is first done in floating point; however if the
c	radius R is found to close to Alpha, or if the "hidden" tests
c	are non conclusive (i.e. test too close to zero), we switch
c	to multiple precision integer arithmetics.
c	The package GMP is used for multiple precision (with a C wrapper)
c
	subroutine alf_tetra(a_xyz,b_xyz,c_xyz,d_xyz,ra,rb,rc,rd,
     1	wa,wb,wc,wd,testr,testa,hull,EPS,SCALE,ALPHA)
c
c	Input:
c		a_xyz,b_xyz	: coordinates of the 4 vertices
c		c_xyz,d_xyz	  that define the tetrahedron
c		wa,wb,wc,wd	: weights of the vertices a,b,c, and d
c				  (i.e. x**2+y**2+z**2-r**2)
c		EPS 		: cutoff value for floating point
c				  filter; if value below, switch
c				  to GMP
c		SCALE		: factor used to convert floating points
c				  to multi precision integers
c		ALPHA		: value of alpha for the alpha shape
c				  (usually 0 in Biology)
c
c	Output:
c		testr_t		: flag for comparison of ALPHA and Rt,
c				  the radius of the circumsphere
c				  of the tetrahedron (.true. if
c				  Rt < Alpha)
c		testr_f		: flag for comparison of ALPHA and
c				  Rf, the radii of the circumspheres
c				  of the four triangles of the tetrahedron
c		testr_e		: flag for comparison of ALPHA and Re,
c				  the radius of the circumpsheres of
c				  the six edges of the tetrahedron
c
c	12/1/2000: this is now obsolete as all flags are stored in a 
c		   single array testr
c
c		testa_f		: test if a vertex is hidden by its
c				  opposite face
c		testa_e		: for each edge, test if a vertex of
c				  the tetrahedron belongs to its circumsphere
c
c	12/1/2000: this is now obsolete as all flags stored in a 
c		   single array testa
c
	integer	i,j,k,ierr
c
	integer	int_testr(15),int_testa(15)
	integer	hull(4)
c
	real*8	Dabc,Dabd,Dacd,Dbcd
	real*8	D1,D2,D3,D4,Det
	real*8	temp1,temp2,temp3,temp4
	real*8	ra,rb,rc,rd,wa,wb,wc,wd
	real*8	r_t,num,den
	real*8  eps,ALPHA,SCALE
	real*8	a_xyz(3),b_xyz(3),c_xyz(3),d_xyz(3)
	real*8	a(4),b(4),c(4),d(4)
	real*8	Sab(3),Sac(3),Sad(3),Sbc(3),Sbd(3),Scd(3)
	real*8	Dab(3),Dac(3),Dad(3),Dbc(3),Dbd(3),Dcd(3)
	real*8	Sa(3),Sb(3),Sc(3),Sd(3)
	real*8	Sam1(3),Sbm1(3),Scm1(3),Sdm1(3)
	real*8	Ta(3),Tb(3),Tc(3),Td(3)
	real*8	Tam1(3),Tbm1(3),Tcm1(3),Tdm1(3)
	real*8	Ua(3),Ub(3),Uc(3),Ud(3)
	real*8	Deter(3)
c
	logical test_a,test_r,second(4)
	logical testr(15),testa(15)
	logical testr_l(15),testa_l(15)
c
c	Initialize local flags
c
	do 100 i = 1,15
		testr_l(i) = testr(i)
		testa_l(i) = testa(i)
100	continue
c
c	first create vectors in 4D space, by adding the weights as the
c	fourth coordinate
c
	do 200 i = 1,3
		a(i) = a_xyz(i)
		b(i) = b_xyz(i)
		c(i) = c_xyz(i)
		d(i) = d_xyz(i)
200	continue
	a(4) = wa
	b(4) = wb
	c(4) = wc
	d(4) = wd
c
c	Perform computation in floating points; if a problem occurs,
c	switch to GMP
c
c	1. Computes all Minors Smn(i+j-2)= M(m,n,i,j) = Det | m(i)  m(j) |
c						            | n(i)  n(j) |
c	for all i in [1,2] and all j in [i+1,3]
c
	do 400 i = 1,2
		do 300 j = i+1,3
			k = i+j-2
			Sab(k) = a(i)*b(j)-a(j)*b(i)
			Sac(k) = a(i)*c(j)-a(j)*c(i)
			Sad(k) = a(i)*d(j)-a(j)*d(i)
			Sbc(k) = b(i)*c(j)-b(j)*c(i)
			Sbd(k) = b(i)*d(j)-b(j)*d(i)
			Scd(k) = c(i)*d(j)-c(j)*d(i)
300		continue
400	continue
c
c	Now compute all Minors 
c		Sq(i+j-2) = M(m,n,p,i,j,0) = Det | m(i) m(j) 1 |
c		       			         | n(i) n(j) 1 |
c						 | p(i) p(j) 1 |
c
c	and all Minors
c		Det(i+j-2) = M(m,n,p,q,i,j,4,0) = Det | m(i) m(j) m(4) 1 |
c						      | n(i) n(j) n(4) 1 |
c						      | p(i) p(j) p(4) 1 |
c						      | q(i) q(j) q(4) 1 |
c
c	m,n,p,q are the four vertices of the tetrahedron, i and j correspond
c	to two of the coordinates of the vertices, and m(4) refers to the
c	"weight" of vertices m
c
	do 500 i = 1,3
		Sa(i) = Scd(i) - Sbd(i) + Sbc(i)
		Sb(i) = Scd(i) - Sad(i) + Sac(i)
		Sc(i) = Sbd(i) - Sad(i) + Sab(i)
		Sd(i) = Sbc(i) - Sac(i) + Sab(i)
		Sam1(i) = -Sa(i)
		Sbm1(i) = -Sb(i)
		Scm1(i) = -Sc(i)
		Sdm1(i) = -Sd(i)
500	continue
c
	do 600 i = 1,3
		Deter(i) = a(4)*Sa(i)-b(4)*Sb(i)+c(4)*Sc(i)-d(4)*Sd(i)
600	continue
c
c	Now compute the determinant needed to compute the radius of the
c	circumsphere of the tetrahedron :
c
c		D1 = Minor(a,b,c,d,4,2,3,0)
c		D2 = Minor(a,b,c,d,1,3,4,0)
c		D3 = Minor(a,b,c,d,1,2,4,0)
c		D4 = Minor(a,b,c,d,1,2,3,0)
c
	D1 = Deter(3)
	D2 = Deter(2)
	D3 = Deter(1)
	D4 = a(1)*Sa(3)-b(1)*Sb(3)+c(1)*Sc(3)-d(1)*Sd(3)
c
c	Now compute all minors:
c		Dmnp = Minor(m,n,p,1,2,3) = Det | m(1) m(2) m(3) |
c						| n(1) n(2) n(3) |
c						| p(1) p(2) p(3) |
c
	Dabc = a(1)*Sbc(3)-b(1)*Sac(3) + c(1)*Sab(3)
	Dabd = a(1)*Sbd(3)-b(1)*Sad(3) + d(1)*Sab(3)
	Dacd = a(1)*Scd(3)-c(1)*Sad(3) + d(1)*Sac(3)
	Dbcd = b(1)*Scd(3)-c(1)*Sbd(3) + d(1)*Sbc(3)
c
c	We also need :
c		Det = Det | m(1) m(2) m(3) m(4) |
c			  | n(1) n(2) n(3) n(4) |
c			  | p(1) p(2) p(3) p(4) |
c			  | q(1) q(2) q(3) q(4) |
c
	Det = -a(4)*Dbcd + b(4)*Dacd -c(4)*Dabd + d(4)*Dabc
c
c	The radius of the circumsphere of the weighted tetrahedron is then:
c
	num = (D1*D1 + D2*D2 + D3*D3 + 4*D4*Det)
	den = (4*D4*D4)
	r_t = num/den
c
c	If this radius is too close to the value of ALPHA, we switch to
c	GMP
c
	if(abs(num-alpha).lt.eps) goto 1200
c
c	The spectrum for a tetrahedron is [R_t Infinity[. If ALPHA is in
c	that interval, the tetrahedron is part of the alpha shape, otherwise
c	it is discarded
c	If tetrahedron is part of the alpha shape, then the 4 triangles,
c	the 6 edges and the four vertices are also part of the alpha
c	complex
c
	if(ALPHA.gt.r_t) then
c
		do 700 i = 1,15
c
			testr(i) = .True.
c
c			Artificially set all triangles and pairs as
c			"attached", so that they will not be checked
c			again
c
			testa(i) = .True.
c
700		continue
c
		return
c
	else
		testr_l(1) = .False.
	endif
c
c	Now check all four faces of the tetrahedra:
c	We look both if the fourth vertex is "hidden" by the face of 
c	interest, and then we compute the radius of the circumsphere
c	of the face
c
c	We first need the minors:
c
c	Tm(i) = Minor(n,p,q,i,4,0) = det | n(i)  n(4)  1 |
c					 | p(i)  p(4)  1 |
c					 | q(i)  q(4)  1 |
c
c	Um(i) = Minor(n,p,q,i,j,4) = det | n(i) n(j) n(4) |
c					 | p(i) p(j) p(4) |
c					 | q(i) q(j) q(4) |
c
	do 800 i = 1,3
		Dab(i) = a(i) - b(i)
		Dac(i) = a(i) - c(i)
		Dad(i) = a(i) - d(i)
		Dbc(i) = b(i) - c(i)
		Dbd(i) = b(i) - d(i)
		Dcd(i) = c(i) - d(i)
800	continue
c
	do 900 i = 1,3
		Ta(i) = -b(4)*Dcd(i)+c(4)*Dbd(i) - d(4)*Dbc(i)
		Tb(i) = -a(4)*Dcd(i)+c(4)*Dad(i) - d(4)*Dac(i)
		Tc(i) = -a(4)*Dbd(i)+b(4)*Dad(i) - d(4)*Dab(i)
		Td(i) = -a(4)*Dbc(i)+b(4)*Dac(i) - c(4)*Dab(i)
		Ua(i) = b(4)*Scd(i) - c(4)*Sbd(i) + d(4)*Sbc(i)
		Ub(i) = a(4)*Scd(i) - c(4)*Sad(i) + d(4)*Sac(i)
		Uc(i) = a(4)*Sbd(i) - b(4)*Sad(i) + d(4)*Sab(i)
		Ud(i) = a(4)*Sbc(i) - b(4)*Sac(i) + c(4)*Sab(i)
		Tam1(i) = -Ta(i)
		Tbm1(i) = -Tb(i)
		Tcm1(i) = -Tc(i)
		Tdm1(i) = -Td(i)
900	continue
c
c	First check face abc:
c
	temp1 = -D1
	temp2 = -D2
	temp3 = -D3
	temp4 = -D4
c
	if(.not.testa_l(2)) then
		test_r = testr_l(2)
		test_a = testa_l(2)
		second(1) = test_r
		call check_triangle(Sd,Td,Ud,D1,D2,D3,D4,Dabc,
     1		test_r,test_a,ALPHA,eps,ierr)
		if(ierr.eq.1) goto 1200
		testa_l(2) = testa_l(2).or.test_a
		testr_l(2) = testr_l(2).or.test_r
		if(testa_l(2)) testr_l(2) = .False.
		if(hull(1).eq.1) second(1)=testr_l(2)
	endif
c
c	Now check face abd
c
	if(.not.testa_l(3)) then
		test_r = testr_l(3)
		test_a = testa_l(3)
		second(2) = test_r
		call check_triangle(Sc,Tc,Uc,temp1,temp2,
     1		temp3,temp4,Dabd,test_r,test_a,ALPHA,eps,ierr)
		if(ierr.eq.1) goto 1200
		testa_l(3) = testa_l(3).or.test_a
		testr_l(3) = testr_l(3).or.test_r
		if(testa_l(3)) testr_l(3) = .False.
		if(hull(2).eq.1) second(2)=testr_l(3)
	endif
c
c	Now check face acd
c
	if(.not.testa_l(4)) then
		test_r = testr_l(4)
		test_a = testa_l(4)
		second(3) = test_r
		call check_triangle(Sb,Tb,Ub,D1,D2,D3,D4,Dacd,
     1		test_r,test_a,ALPHA,eps,ierr)
		if(ierr.eq.1) goto 1200
		testa_l(4) = testa_l(4).or.test_a
		testr_l(4) = testr_l(4).or.test_r
		if(testa_l(4)) testr_l(4) = .False.
		if(hull(3).eq.1) second(3)=testr_l(4)
	endif
c
c	Now check face bcd
c
	if(.not.testa_l(5)) then
		test_r = testr_l(5)
		test_a = testa_l(5)
		second(4) = test_r
		call check_triangle(Sa,Ta,Ua,temp1,temp2,temp3,
     1		temp4,Dbcd,test_r,test_a,ALPHA,eps,ierr)
		if(ierr.eq.1) goto 1200
		testa_l(5) = testa_l(5).or.test_a
		testr_l(5) = testr_l(5).or.test_r
		if(testa_l(5)) testr_l(5) = .False.
		if(hull(4).eq.1) second(4)=testr_l(5)
	endif
c
c
c	Now consider edges
c
c	Start by checking each triangle: if it belongs to the alpha shape,
c	so does its 3 edges:
c
	if(testr_l(2).and.second(1)) then
		testr_l(6) = .True.
		testr_l(7) = .True.
		testr_l(9) = .True.
		testa_l(6) = .True.
		testa_l(7) = .True.
		testa_l(9) = .True.
	endif
c
	if(testr_l(3).and.second(2)) then
		testr_l(6) = .True.
		testr_l(8) = .True.
		testr_l(10) = .True.
		testa_l(6) = .True.
		testa_l(8) = .True.
		testa_l(10) = .True.
	endif
c
	if(testr_l(4).and.second(3)) then
		testr_l(7) = .True.
		testr_l(8) = .True.
		testr_l(11) = .True.
		testa_l(7) = .True.
		testa_l(8) = .True.
		testa_l(11) = .True.
	endif
c
	if(testr_l(5).and.second(4)) then
		testr_l(9) = .True.
		testr_l(10) = .True.
		testr_l(11) = .True.
		testa_l(9) = .True.
		testa_l(10) = .True.
		testa_l(11) = .True.
	endif
c
c	Now check each edges
c
	if(.not.testa_l(6)) then
		test_r = testr_l(6)
		test_a = testa_l(6)
		call check_edge(a_xyz,b_xyz,wa,wb,Dab,
     1		Sab,Sd,Sc,Td,Tc,test_r,test_a,alpha,eps,ierr)
		if(ierr.eq.1) goto 1200
		testa_l(6) = testa_l(6).or.test_a
		testr_l(6) = testr_l(6).or.test_r
		if(testa_l(6)) testr_l(6) = .False.
	endif
c
	if(.not.testa_l(7)) then
		test_r = testr_l(7)
		test_a = testa_l(7)
		call check_edge(a_xyz,c_xyz,wa,wc,Dac,Sac,
     1		Sdm1,Sb,Tdm1,Tb,test_r,test_a,alpha,eps,ierr)
		if(ierr.eq.1) goto 1200
		testa_l(7) = testa_l(7).or.test_a
		testr_l(7) = testr_l(7).or.test_r
		if(testa_l(7)) testr_l(7) = .False.
	endif
c
	if(.not.testa_l(8)) then
		test_r = testr_l(8)
		test_a = testa_l(8)
		call check_edge(a_xyz,d_xyz,wa,wd,Dad,Sad,
     1		Scm1,Sbm1,Tcm1,Tbm1,test_r,test_a,alpha,eps,ierr)
		if(ierr.eq.1) goto 1200
		testa_l(8) = testa_l(8).or.test_a
		testr_l(8) = testr_l(8).or.test_r
		if(testa_l(8)) testr_l(8) = .False.
	endif
c
	if(.not.testa_l(9)) then
		test_r = testr_l(9)
		test_a = testa_l(9)
		call check_edge(b_xyz,c_xyz,wb,wc,Dbc,Sbc,
     1		Sd,Sa,Td,Ta,test_r,test_a,alpha,eps,ierr)
		if(ierr.eq.1) goto 1200
		testa_l(9) = testa_l(9).or.test_a
		testr_l(9) = testr_l(9).or.test_r
		if(testa_l(9)) testr_l(9) = .False.
	endif
c
	if(.not.testa_l(10)) then
		test_r = testr_l(10)
		test_a = testa_l(10)
		call check_edge(b_xyz,d_xyz,wb,wd,Dbd,Sbd,Sc,Sam1,
     1		Tc,Tam1,test_r,test_a,alpha,eps,ierr)
		if(ierr.eq.1) goto 1200
		testa_l(10) = testa_l(10).or.test_a
		testr_l(10) = testr_l(10).or.test_r
		if(testa_l(10)) testr_l(10) = .False.
	endif
c
	if(.not.testa_l(11)) then
		test_r = testr_l(11)
		test_a = testa_l(11)
		call check_edge(c_xyz,d_xyz,wc,wd,Dcd,Scd,Sb,
     1		Sa,Tb,Ta,test_r,test_a,alpha,eps,ierr)
		if(ierr.eq.1) goto 1200
		testa_l(11) = testa_l(11).or.test_a
		testr_l(11) = testr_l(11).or.test_r
		if(testa_l(11)) testr_l(11) = .False.
	endif
c
c	Set all vertices in the alpha complex
c
	do 1000 i = 12,15
		testr_l(i) = .True.
1000	continue
c
c	Copy back flags
c
	do 1100 i = 1,15
		testr(i) = testr_l(i)
		testa(i) = testa_l(i)
1100	continue
c
	return
c
c	If a precision problem was accountered, the whole procedure
c	is repeated using multiple precision arithmetics, using
c	the package GMP (Gnu Multi Precision), in C:
c
1200	continue
c
c	Get initial flags
c
	do 1300 i = 1,15
		if(testr(i)) then
			int_testr(i) = 1
		else
			int_testr(i) = 0
		endif
		if(testa(i)) then
			int_testa(i) = 1
		else
			int_testa(i) = 0
		endif
1300	continue
c
	call alf_tetra_gmp(a,b,c,d,ra,rb,rc,rd,int_testr,
     1		int_testa,SCALE,ALPHA)
c
	do 1400 i = 1,15
		if(int_testr(i).eq.1) then
			testr(i) = .True.
		else
			testr(i) = .False.
		endif
		if(int_testa(i).eq.1) then
			testa(i) = .True.
		else
			testa(i) = .False.
		endif
1400	continue
c
	return
	end
c
c	Check_triangle.f	Version 1 11/22/2000	Patrice Koehl
c
	subroutine check_triangle(S,T,U,Det1,Det2,Det3,Deter,De3,
     1	testr,testa,alpha,eps,ierr)
c
c	Input:
c
c	For the three points a,b,c that form the triangles, the program
c	needs as input the following determinants:
c
c	S(i+j-2) = Minor(a,b,c,i,j,0)= det | a(i)  a(j)  1 |
c					   | b(i)  b(j)  1 |
c					   | c(i)  c(j)  1 |
c
c	T(i) = Minor(a,b,c,i,4,0) = det | a(i)  a(4)  1 |
c					| b(i)  b(4)  1 |
c					| c(i)  c(4)  1 |
c
c	U(i) = Minor(a,b,c,i,j,4) = det | a(i) a(j) a(4) |
c					| b(i) b(j) b(4) |
c					| c(i) c(j) c(4) |
c
c
c	Det1 = Minor(a,b,c,d,2,3,4,0)
c	Det2 = Minor(a,b,c,d,1,3,4,0)
c	Det3 = Minor(a,b,c,d,1,2,4,0)
c	Deter= Minor(a,b,c,d,1,2,3,0)
c	De3  = Minor(a,b,c,1,2,3)
c
c	Output:
c
c	testr	: flag set to 1 if ALPHA is larger than rho, the radius
c		  of the circumsphere of the triangle
c	testa	: flag set to 1 if the fourth point d of the tetrahedron
c		  that contains the triangle {a,b,c} is inside the
c		  circumsphere of {a,b,c}
c
c	The program tests for problem with floating points (i.e. some
c	results smaller than EPS, in which case the sign of the
c	expression cannot be defined). If problems, IERR returns as 1,
c	and the program will then switch to LIA
c
	integer	ierr
c
	real*8	Det1,Det2,Det3,Deter,De3
	real*8	d0,d1,d2,d3,d4,test
	real*8	eps,alpha
	real*8	sums2,num
	real*8	S(3),T(3),U(3)
c
	logical testr,testa
c
	ierr = 0
c
	sums2 = S(1)*S(1)+S(2)*S(2)+S(3)*S(3)
c
c	First check if the face is "attached" to the four vertex of the
c	parent tetrahedron
c	
c	Skip test if triangle was already found attached
c
	if(testa) goto 100
c
	test =  sums2 *(Det1*S(3)+Det2*S(2)+Det3*S(1)-2*Deter*De3) 
c
c	Check for problems, in which case should be LIA
c
	if(abs(test).lt.eps) then
		ierr = 1
		return
	endif
c
c	If no problem, set testa to true to t > 0
c
	if(test.gt.0) then
		testa = .TRUE.
	else
		testa = .FALSE.
	endif
c
100	continue
c
c	Now compute the radius of the smallest circumsphere of the triangle
c	{a,b,c}
c
	if(testr) return
c
	d0 = 4*sums2
c
	d1 = -2*(S(2)*T(3)+S(1)*T(2) -2*De3*S(3))
	d2 =  2*(S(1)*T(1)-S(3)*T(3) -2*De3*S(2))
	d3 =  2*(S(3)*T(2)+S(2)*T(1) +2*De3*S(1))
	d4 = -4*(S(1)*U(1)+S(2)*U(2) + S(3)*U(3) - 2*De3*De3)
c
	num  = (d1*d1+d2*d2+d3*d3 - d0*d4)
c
	if(abs(alpha-num).lt.eps) then
		ierr = 1
		return
	endif
c
	if(alpha.gt.num) then
		testr = .TRUE.
	else
		testr = .FALSE.
	endif
c
	return
	end
c
c	check_edge.f		Version 1 11/25/2000	Patrice Koehl
c
	subroutine check_edge(a,b,wa,wb,Dab,Sab,Sc,Sd,Tc,Td,
     1	testr,testa,alpha,eps,ierr)
c
c	This subroutine check if an edge of a tetrahedron is "attached"
c	to one of the two remaining vertices of the tetrahedron (i.e.
c	if one of the two vertices belong to the smallest circumsphere
c	of the edge). For that, it needs:
c
c	Input:
c		a,b	: coordinate of the two vertices defining the edge
c		wa,wb	: weights of these two vertices
c		Dab	: minor(a,b,i,0) for all i=1,2,3
c		Sab	: minor(a,b,i,j) for i = 1,2 and j =i+1,3
c		Sc	: minor(a,b,c,i,j,0) for i=1,2 and j = i+1,3
c		Sd	: minor(a,b,d,i,j,0) for i=1,2 and j = i+1,3
c			  c and d are the two other vertices of the
c			  tetrahedron
c		Tc	: minor(a,b,c,i,4,0) for i = 1,2,3
c		Td	: minor(a,b,d,i,4,0) for i = 1,2,3
c		alpha	: value of alpha considered
c		eps	: precision: if a floating point test lead to a
c			  value below this precision, computation
c			  switches to LIA
c	Ouput:
c		testr	: flag that defines if edge belongs to the
c			  Alpha complex
c		testa	: flag that defines if edge is attached or not
c		ierr	: flag to decide if  computation must be done
c			  in LIA
c
	integer	i,ierr
c
	real*8	wa,wb,d0,d1,d2,d3,d4,d5
	real*8	alpha,eps,num,den,rho2,dtest
	real*8	r_11,r_22,r_33,r_14,r_313,r_212,diff
	real*8	a(3),b(3)
	real*8	Sab(3),Dab(3),Sc(3),Sd(3),Tc(3),Td(3)
	real*8	res(0:3,1:4),res2_c(3,4),res2_d(3,4)
c
	logical testr,testa
c
	ierr = 0
c
c	Formulas have been derived by projection on 4D space,
c	which requires some precaution when some coordinates are
c	equal.
c
	res(0,4) = wa-wb
c
	if(a(1).ne.b(1)) then
		do 100 i = 1,3
			res(0,i) = Dab(i)
			res2_c(i,4) = Tc(i)
			res2_d(i,4) = Td(i)
			res(i,4) = a(i)*wb-b(i)*wa
100		continue
		res(1,2) = Sab(1)
		res(1,3) = Sab(2)
		res(2,3) = Sab(3)
		res2_c(1,2) = Sc(1)
		res2_c(1,3) = Sc(2)
		res2_c(2,3) = Sc(3)
		res2_d(1,2) = Sd(1)
		res2_d(1,3) = Sd(2)
		res2_d(2,3) = Sd(3)
	elseif(a(2).ne.b(2)) then
		res(0,1) = Dab(2)
		res(0,2) = Dab(3)
		res(0,3) = Dab(1)
		res(1,2) = Sab(3)
		res(1,3) = -Sab(1)
		res(2,3) = -Sab(2)
		res(1,4) = a(2)*wb-b(2)*wa
		res(2,4) = a(3)*wb-b(3)*wa
		res(3,4) = a(1)*wb-b(1)*wa
		res2_c(1,2) = Sc(3)
		res2_c(1,3) = -Sc(1)
		res2_c(2,3) = -Sc(2)
		res2_d(1,2) = Sd(3)
		res2_d(1,3) = -Sd(1)
		res2_d(2,3) = -Sd(2)
		res2_c(1,4) = Tc(2)
		res2_c(2,4) = Tc(3)
		res2_c(3,4) = Tc(1)
		res2_d(1,4) = Td(2)
		res2_d(2,4) = Td(3)
		res2_d(3,4) = Td(1)
	elseif(a(3).ne.b(3)) then
		res(0,1) = Dab(3)
		res(0,2) = Dab(1)
		res(0,3) = Dab(2)
		res(1,2) = -Sab(2)
		res(1,3) = -Sab(3)
		res(2,3) = Sab(1)
		res(1,4) = a(3)*wb-b(3)*wa
		res(2,4) = a(1)*wb-b(1)*wa
		res(3,4) = a(2)*wb-b(2)*wa
		res2_c(1,2) = -Sc(2)
		res2_c(1,3) = -Sc(3)
		res2_c(2,3) = Sc(1)
		res2_d(1,2) = -Sd(2)
		res2_d(1,3) = -Sd(3)
		res2_d(2,3) = Sd(1)
		res2_c(1,4) = Tc(3)
		res2_c(2,4) = Tc(1)
		res2_c(3,4) = Tc(2)
		res2_d(1,4) = Td(3)
		res2_d(2,4) = Td(1)
		res2_d(3,4) = Td(2)
	else
		write(6,*) 'Problem in hidden1: edges defined from',
     1		' a single point'
		stop
	endif
c
	r_11 = res(0,1)*res(0,1)
	r_22 = res(0,2)*res(0,2)
	r_33 = res(0,3)*res(0,3)
	r_14 = res(0,1)*res(0,4)
	r_313 = res(0,3)*res(1,3)
	r_212 = res(0,2)*res(1,2)
	diff = res(0,3)*res(1,2)-res(0,2)*res(1,3)
c
c	First compute radius of circumsphere
c
	d0 = -2*res(0,1)*(r_11+r_22+r_33)
	d1 = res(0,1)*(2*(r_313 + r_212)-r_14)
	d2 = -2*res(1,2)*(r_11+r_33) - res(0,2)*(r_14-2*r_313)
	d3 = -2*res(1,3)*(r_11+r_22) -res(0,3)*(r_14-2*r_212)
	d4 = 2*res(0,1)*(res(0,1)*res(1,4)+res(0,2)*res(2,4)+
     1		res(0,3)*res(3,4)) +4*(res(2,3)*diff
     2		    - res(0,1)*(res(1,2)*res(1,2)+res(1,3)*res(1,3)))
c
	num = d1*d1 + d2*d2 + d3*d3 - d0*d4
	den = d0*d0
c
	rho2 = num/den
c
	if(abs(alpha-rho2).lt.eps) then
		ierr = 1
		return
	endif
c
	if(alpha.gt.rho2) then
		testr = .TRUE.
	else
		testr = .FALSE.
	endif
c
c	Now check attachement: first test vertex C
c
c	skip if already found attached
c
	if(testa) goto 200
c
	d5 = res(0,1)*(res(0,1)*res2_c(1,4)+res(0,2)*res2_c(2,4)
     1		+res(0,3)*res2_c(3,4)-2*(res(1,3)*res2_c(1,3)
     2		+res(1,2)*res2_c(1,2))) + 2*res2_c(2,3)*diff
c
	dtest = d0*d5
c
	if(abs(dtest).lt.eps) then
		ierr = 1
		return
	endif
c
c	If no problem, set testa to true if t < 0
c
	if(dtest.lt.0) then
		testa = .TRUE.
	else
		testa = .FALSE.
	endif
c
200	continue
c
	if(testa) return
c
	d5 = res(0,1)*(res(0,1)*res2_d(1,4)+res(0,2)*res2_d(2,4)
     1		+res(0,3)*res2_d(3,4)-2*(res(1,3)*res2_d(1,3)
     2		+res(1,2)*res2_d(1,2))) + 2*res2_d(2,3)*diff
c
	dtest = d0*d5
c
	if(abs(dtest).lt.eps) then
		ierr = 1
		return
	endif
c
c	If no problem, set testa to true if t < 0
c
	if(dtest.lt.0) then
		testa = .TRUE.
	else
		testa = .FALSE.
	endif
c
	return
	end
