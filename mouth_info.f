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
c	Find_mouth.f	Version 1 4/10/2002	Patrice Koehl
c
c	This subroutine defines all mouths of the pockets of the
c	protein considered
c
	subroutine find_mouth
c
	include 'pocket.h'
c
	integer	i,j,k,l,n,itrig,itest
	integer ntetra,ntrig,nmouth,npocket
	integer	iset1,iset2
	integer	pair1_i,pair2_i,pair3_i
	integer	pair1,pair2,pair3
c
	integer	trig(3,ntrig_max)
c
	integer*1 trig_status(ntrig_max)
	integer*1 trig_coef(ntrig_max)
	integer trig_link(3,ntrig_max)
c
	integer	tetra(4,ntetra_max),tetra_neighbour(4,ntetra_max)
	integer tetra_link(4,ntetra_max)
	integer	tetra_pocket(ntetra_max)
c
	integer trig_root(ntrig_max),trig_pointer(ntrig_max)
	integer	mouth_root(nmouth_max),mouth_count(nmouth_max)
	integer	mouth_pocket(nmouth_max)
c
	common	/trig_zone/	ntrig,trig
	common  /trig_stat/	trig_status,trig_coef
	common  /tetra_zone/	ntetra,tetra,tetra_neighbour
	common  /links/		tetra_link,trig_link
	common  /pockets/	npocket,tetra_pocket
	common  /mouths/	nmouth,mouth_pocket,trig_root
c
	do 50 i = 1,ntrig
		trig_root(i) = 0
50	continue
c
	nmouth = 0
c
	do 200 i = 1,ntetra
		if(tetra_pocket(i).eq.0) goto 200
		do 100 j = 1,4
			n = tetra_neighbour(j,i)
			itrig=tetra_link(j,i)
			if(n.eq.0) then
				if(trig_status(itrig).eq.0) then
					nmouth= nmouth + 1
					trig_root(itrig)=nmouth
					trig_pointer(itrig)=0
					mouth_root(nmouth)=itrig
					mouth_count(nmouth)=1
					mouth_pocket(nmouth)=
     1					tetra_pocket(i)
				endif
			else
				if(tetra_pocket(n).ne.tetra_pocket(i).
     1				and.trig_status(itrig).eq.0) then
					nmouth= nmouth + 1
					trig_root(itrig)=nmouth
					trig_pointer(itrig)=0
					mouth_root(nmouth)=itrig
					mouth_count(nmouth)=1
					mouth_pocket(nmouth)=
     1					tetra_pocket(i)
				endif
			endif
100		continue
200	continue
c
	do 700 i = 1,ntrig
c
		iset1 = trig_root(i)
		if(iset1.eq.0) goto 700
		pair1_i = trig_link(1,i)
		pair2_i = trig_link(2,i)
		pair3_i = trig_link(3,i)
c
		do 600 j = i+1,ntrig
			iset2 = trig_root(j)
			if(iset2.eq.0) goto 600
			iset1 = trig_root(i)
			if(iset2.eq.iset1) goto 600
			pair1 = trig_link(1,j)
			pair2 = trig_link(2,j)
			pair3 = trig_link(3,j)
			itest = 0
			if(pair1_i.eq.pair1) then
				itest = 1
				goto 300
			endif
			if(pair1_i.eq.pair2) then
				itest = 1
				goto 300
			endif
			if(pair1_i.eq.pair3) then
				itest = 1
				goto 300
			endif
			if(pair2_i.eq.pair1) then
				itest = 1
				goto 300
			endif
			if(pair2_i.eq.pair2) then
				itest = 1
				goto 300
			endif
			if(pair2_i.eq.pair3) then
				itest = 1
				goto 300
			endif
			if(pair3_i.eq.pair1) then
				itest = 1
				goto 300
			endif
			if(pair3_i.eq.pair2) then
				itest = 1
				goto 300
			endif
			if(pair3_i.eq.pair3) itest = 1
300			continue
			if(itest.eq.1) then
			   if(mouth_count(iset2).lt.
     1			   mouth_count(iset1)) then
				l = mouth_root(iset2)
				do 400 k = 1,mouth_count(iset2)-1
					trig_root(l) = iset1
					l = trig_pointer(l)
400				continue
				trig_root(l) = iset1
				trig_pointer(l)= mouth_root(iset1)
				mouth_root(iset1) = mouth_root(iset2)
				mouth_count(iset1)=mouth_count(iset1)+
     1				mouth_count(iset2)
				mouth_count(iset2) = 0
			   else
				l = mouth_root(iset1)
				do 500 k = 1,mouth_count(iset1)-1
					trig_root(l) = iset2
					l = trig_pointer(l)
500				continue
				trig_root(l) = iset2
				trig_pointer(l)= mouth_root(iset2)
				mouth_root(iset2) = mouth_root(iset1)
				mouth_count(iset2)=mouth_count(iset2)+
     1				mouth_count(iset1)
				mouth_count(iset1) = 0
			   endif
			endif
600		continue
c
700	continue
c
	do 800 i = 1,ntrig
		trig_root(i) = 0
800	continue
c
	k = 0
	do 1000 i = 1,nmouth
		if(mouth_count(i).eq.0) goto 1000
		k = k + 1
		mouth_pocket(k)=mouth_pocket(i)
		l = mouth_root(i)
		do 900 j = 1,mouth_count(i)
			trig_root(l) = k
			l = trig_pointer(l)
900		continue
1000	continue
c
	nmouth = k
c
	return
	end
c
c	Measure_mouth.f		Version 1 8/4/1998	Patrice Koehl
c
c	Based on the alpha shape, this subroutine computes the surface
c	of the mouths of each pocket (cavity) in the space filling
c	diagram
c
	subroutine measure_mouth(surf_mouth,mouth_nvertex,mouth_vertex)
c
	include 'pocket.h'
c
c	Misc. variables
c
	integer	i,j,k,m
	integer	imouth,idx
	integer ntrig
	integer npoints,nvertex
c
c	Information on mouth
c
	integer	nmouth,mouth_nvertex
	integer	trig_mouth(ntrig_max)
	integer	mouth_pocket(nmouth_max)
	integer	mouth_vertex(npointmax)
	integer	vertex(npointmax)
c
c	Information on triangles
c
	integer trig(3,ntrig_max)
c
c	Information on vertices
c
	integer redinfo(npointmax)
c
	real*8	ra,rb,rc
	real*8	surfa,surf,s
	real*8	ang1,ang2,ang3
	real*8	dist1,dist2,dist3
	real*8	a(3),b(3),c(3)
	real*8	surf_mouth(nmouth_max)
c
	real*8  radius(npointmax)
	real*8	weight(npointmax)
	real*8  coord(3*npointmax)
c
	common	/xyz_vertex/	coord,radius,weight
	common	/trig_zone/	ntrig,trig
	common	/vertex_zone/	npoints,nvertex,redinfo
	common /mouths/		nmouth,mouth_pocket,trig_mouth
c	
	include 'constants.h'
c
	mouth_nvertex = 0
c
	do 2000 imouth = 1,nmouth
c
	    surf_mouth(imouth) = 0
c
	    mouth_nvertex = mouth_nvertex + 1
	    mouth_vertex(mouth_nvertex) = - imouth
c
	    do 200 i = 1,npoints
		vertex(i) = 0
200	    continue
c
	    do 1700 idx = 1,ntrig
c
		if(trig_mouth(idx).ne.imouth) goto 1700
c
		i = trig(1,idx)
		j = trig(2,idx)
		k = trig(3,idx)
c
		ra = radius(i)
		rb = radius(j)
		rc = radius(k)
c
		do 400 m = 1,3
			a(m) = coord(3*(i-1)+m)
			b(m) = coord(3*(j-1)+m)
			c(m) = coord(3*(k-1)+m)
400		continue
c
		call triangle_ang(a,b,c,ang1,ang2,ang3,dist1,dist2,
     1		dist3)
c
		s = (dist1+dist2+dist3)/2
		surf = s*(s-dist1)*(s-dist2)*(s-dist3)
		surf = sqrt(surf)
		surf_mouth(imouth) = surf_mouth(imouth) + surf
c
c	First check all vertices of the tetrahedron
c
		if(redinfo(i).eq.0) then
			vertex(i) = imouth
			surfa = pi*ra*ra*ang1
			surf_mouth(imouth)=surf_mouth(imouth)-surfa
		endif
c
		if(redinfo(j).eq.0) then
			vertex(j) = imouth
			surfa = pi*rb*rb*ang2
			surf_mouth(imouth)=surf_mouth(imouth)-surfa
		endif
c
		if(redinfo(k).eq.0) then
			vertex(k) = imouth
			surfa = pi*rc*rc*ang3
			surf_mouth(imouth)=surf_mouth(imouth)-surfa
		endif
c
c	Now check all edges
c
		if(dist1.le.ra+rb) then
			call twocircle(ra,rb,dist1,surfa)
			surf_mouth(imouth)=surf_mouth(imouth) +
     1			0.5*surfa
		endif
c
		if(dist2.le.ra+rc) then
			call twocircle(ra,rc,dist2,surfa)
			surf_mouth(imouth)=surf_mouth(imouth) +
     1			0.5*surfa
		endif
c
		if(dist3.le.rb+rc) then
			call twocircle(rb,rc,dist3,surfa)
			surf_mouth(imouth)=surf_mouth(imouth) +
     1			0.5*surfa
		endif
c
		if(surf_mouth(imouth).lt.0) surf_mouth(imouth) = 0
c
1700	   continue
c
	   do 1800 i = 1,npoints
		if(vertex(i).ne.0) then
			mouth_nvertex=mouth_nvertex+1
			mouth_vertex(mouth_nvertex)=i
		endif
		vertex(i) = 0
1800	   continue
c
2000	continue
c
	return
	end
c
c	triangle_ang.f		Version 1 5/20/2002	Patrice Koehl
c
c	This subroutine computes the three angles of a triangle (ABC)
c
	subroutine triangle_ang(a,b,c,ang1,ang2,ang3,dist1,dist2,dist3)
c
	real*8	x,val,eps
	real*8	ang1,ang2,ang3
	real*8	dist1,dist2,dist3
	real*8	pi,twopi
c
	real*8	a(3),b(3),c(3)
	real*8	u1(3),u2(3),u3(3)
c
	eps = 1.e-8
	pi = acos(-1.)
	twopi = pi*2.
c
	call diffvect(a,b,u1)
	call normvect(u1,dist1)
c
	call diffvect(a,c,u2)
	call normvect(u2,dist2)
c
	call diffvect(b,c,u3)
	call normvect(u3,dist3)
c
	call dotvect(u1,u2,val)
	x = val/(dist1*dist2)
	if(x.ge.1.) x = 1.-eps
	if(x.le.-1) x = -1 + eps
	ang1 = acos(x)/twopi
c
	call dotvect(u1,u3,val)
	x = -val/(dist1*dist3)
	if(x.ge.1.) x = 1.-eps
	if(x.le.-1) x = -1 + eps
	ang2 = acos(x)/twopi
c
	ang3 = 0.5-ang1-ang2
c
	return
	end
c
c	twocircle.f	Version 1 5/20/2002	Patrice Koehl
c
c	This subroutine computes the surface area of the intersection of
c	two disks
c
	subroutine twocircle(ra,rb,dist,surf)
c
	real*8	ra,ra2,rb,rb2,dist,dist2,surf
	real*8	val,eps,ang1,ang2,s
c
	eps = 1.e-8
c
	ra2 = ra*ra
	rb2 = rb*rb
	dist2 = dist*dist
c
	val = (dist2+ra2-rb2)/(2*dist*ra)
	if(val.ge.1.) val = 1.-eps
	if(val.lt.-1.) val = -1+eps
	ang1 = acos(val)
c
	val = (dist2+rb2-ra2)/(2*dist*rb)
	if(val.ge.1.) val = 1.-eps
	if(val.lt.-1.) val = -1+eps
	ang2 = acos(val)
c
	s = (-dist+ra+rb)*(dist+ra-rb)*(dist-ra+rb)*(dist+ra+rb)
	s = 0.5*sqrt(s)
c
	surf = ra2*ang1 + rb2*ang2 - s
c
	return
	end
