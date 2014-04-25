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
c	flow.f	Version 1 11/24/2000	Patrice Koehl
c
c	This subroutine computes the flow coming from a tetrahedron
c
c	To check the flow through a face (A,B,C) of a tetrahedron ABCD,
c	we check the position of the orthocenter O of ABCD with respect
c	to the face (A,B,C): if O and D (the "other" vertex of the
c	tetrahedron) are on the same side, the flow is "in", and "out"
c	otherwise. This test is done by checking if D is "hidden" by
c	(A,B,C). A similar test is used in alpha_fixed
c
c	Computation extends to all four faces of the tetrehedron
c
c	Computation is first done in floating point; however if the "hidden" 
c	tests are non conclusive (i.e. test too close to zero), we switch
c	to multiple precision integer arithmetics.
c	The package GMP is used for multiple precision (with a C wrapper)
c
	subroutine flow(a_xyz,b_xyz,c_xyz,d_xyz,ra,rb,rc,rd,
     1	wa,wb,wc,wd,testa,EPS,SCALE)
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
c
c	Output:
c
c		testa(4)	: test if a vertex is hidden by its
c				  opposite face
c
	integer	i,j,k,ierr
c
	integer	itesta,testa(4)
c
	real*8	Dabc,Dabd,Dacd,Dbcd
	real*8	D1,D2,D3,D4
	real*8	temp1,temp2,temp3,temp4
	real*8	ra,rb,rc,rd,wa,wb,wc,wd
	real*8  eps,SCALE
	real*8	a_xyz(3),b_xyz(3),c_xyz(3),d_xyz(3)
	real*8	a(4),b(4),c(4),d(4)
	real*8	Sab(3),Sac(3),Sad(3),Sbc(3),Sbd(3),Scd(3)
	real*8	Sa(3),Sb(3),Sc(3),Sd(3)
	real*8	Deter(3)
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
c	Now check all four faces of the tetrahedra:
c	We look both if the fourth vertex is "hidden" by the face of 
c	interest
c
c	First check face abc:
c
	temp1 = -D1
	temp2 = -D2
	temp3 = -D3
	temp4 = -D4
c
	if(testa(4).eq.0) then
		call check_facet(Sd,D1,D2,D3,D4,Dabc,
     1		itesta,eps,ierr)
		if(ierr.eq.1) goto 1200
		testa(4) = itesta
	endif
c
c	Now check face abd
c
	if(testa(3).eq.0) then
		call check_facet(Sc,temp1,temp2,
     1		temp3,temp4,Dabd,itesta,eps,ierr)
		if(ierr.eq.1) goto 1200
		testa(3) = itesta
	endif
c
c	Now check face acd
c
	if(testa(2).eq.0) then
		call check_facet(Sb,D1,D2,D3,D4,Dacd,
     1		itesta,eps,ierr)
		if(ierr.eq.1) goto 1200
		testa(2) = itesta
	endif
c
c	Now check face bcd
c
	if(testa(1).eq.0) then
		call check_facet(Sa,temp1,temp2,temp3,
     1		temp4,Dbcd,itesta,eps,ierr)
		if(ierr.eq.1) goto 1200
		testa(1) = itesta
	endif
c
	return
c
c	If a precision problem was accountered, the whole procedure
c	is repeated using multiple precision arithmetics, using
c	the package GMP (Gnu Multi Precision), in C:
c
1200	continue
c
	call tetra_flow_gmp(a,b,c,d,ra,rb,rc,rd,testa,SCALE)
c
	return
	end
c
c	Check_facet.f	Version 1 11/22/2000	Patrice Koehl
c
	subroutine check_facet(S,Det1,Det2,Det3,Deter,De3,
     1	testa,eps,ierr)
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
c	Det1 = Minor(a,b,c,d,2,3,4,0)
c	Det2 = Minor(a,b,c,d,1,3,4,0)
c	Det3 = Minor(a,b,c,d,1,2,4,0)
c	Deter= Minor(a,b,c,d,1,2,3,0)
c	De3  = Minor(a,b,c,1,2,3)
c
c	Output:
c
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
	integer	testa
c
	real*8	Det1,Det2,Det3,Deter,De3
	real*8	test
	real*8	eps
	real*8	sums2
	real*8	S(3)
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
		testa = 1
	else
		testa = 0
	endif
c
	return
	end
