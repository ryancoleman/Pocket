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
c	find_pockets.f	Version 1 4/10/2002	Patrice Koehl
c
c	This subroutine defines all pockets in the protein considered
c
	subroutine find_pockets
c
	include 'pocket.h'
c
	integer	i,j,k,l,m,n,p
	integer npockets,itrig
	integer ntetra
	integer	iset1,iset2
	integer	ndepth,nlist,ilist,nset
c
	integer	depth_list(ntetra_max),depth_to_list(ntetra_max)
	integer	nlist_count(nlistmax),list(ntet_listmax,nlistmax)
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
	integer*1 trig_status(ntrig_max)
	integer*1 trig_coef(ntrig_max)
	integer trig_link(3,ntrig_max)
c
	integer	tetra(4,ntetra_max),tetra_neighbour(4,ntetra_max)
	integer tetra_link(4,ntetra_max)
	integer	tetra_depth(ntetra_max)
c
	integer tetra_root(ntetra_max),tetra_pointer(ntetra_max)
	integer	pocket_root(ntetra_max),pocket_count(ntetra_max)
c
	common  /trig_stat/	trig_status,trig_coef
	common  /tetra_zone/	ntetra,tetra,tetra_neighbour
	common  /tetra_stat/	tetra_status,tetra_orient,
     1				tetra_nindex
	common  /tet_depth/	tetra_depth
	common  /links/		tetra_link,trig_link
	common  /pockets/	npockets,tetra_root
c
	ndepth = 0
	do 100 i = 1,ntetra
		if(tetra_status(i).ne.-1) goto 100
		ndepth = ndepth + 1
		depth_list(ndepth) = tetra_depth(i)
		tetra_root(i) = 0
100	continue
c
	call hpsort_int(depth_list,ndepth)
c
	nlist = 1
	do 200 i = 2,ndepth
		if(depth_list(i).ne.depth_list(nlist)) then
			nlist = nlist + 1
			depth_list(nlist) = depth_list(i)
		endif
200	continue
c
	do 300 i = 1,ntetra
		depth_to_list(i) = 0
300	continue
c
	do 400 i = 1,nlist
		depth_to_list(depth_list(i)) = i
		nlist_count(i) = 0
400	continue
c
	npockets = 0
c
	do 1000 i = 1,ntetra
c
		if(tetra_status(i).ne.-1) goto 1000
c
		ndepth = tetra_depth(i)
		if(ndepth.eq.ntetra+1) goto 1000
c
		ilist = depth_to_list(ndepth)
		nlist_count(ilist) = nlist_count(ilist) + 1
		list(nlist_count(ilist),ilist) = i
c
		nset = depth_to_list(i)
		if(nset.eq.0) goto 1000
c
		do 900 j = 1,nlist_count(nset)
c
			k = list(j,nset)
c
			npockets = npockets + 1
			tetra_root(k) = npockets
			tetra_pointer(k) = 0
			pocket_root(npockets) = k
			pocket_count(npockets) = 1
c
			do 800 l = 1,4
c
				m = tetra_neighbour(l,k)
				if(m.eq.0) goto 800
				if(tetra_status(m).ne.-1) goto 800
c				if(tetra_depth(m).gt.i) goto 800
				if(tetra_root(m).eq.0) goto 800
c
				itrig = tetra_link(l,k)
				if(trig_status(itrig).eq.1) goto 800
c
				iset1 = tetra_root(k)
				iset2 = tetra_root(m)
				if(iset1.eq.iset2) goto 800
c
				if(pocket_count(iset1).gt.
     1				pocket_count(iset2)) then
c
				   n = pocket_root(iset2)
				   do 500 p=1,pocket_count(iset2)-1
					tetra_root(n) = iset1
					n = tetra_pointer(n)
500				   continue
				   tetra_root(n) = iset1
				   tetra_pointer(n) = 
     1					pocket_root(iset1)
				   pocket_root(iset1)=
     1					pocket_root(iset2)
				   pocket_count(iset1)=
     1				   pocket_count(iset1)+
     2					pocket_count(iset2)
				   pocket_count(iset2) = 0
c
				else
c
				   n = pocket_root(iset1)
				   do 600 p=1,pocket_count(iset1)-1
					tetra_root(n) = iset2
					n = tetra_pointer(n)
600				   continue
				   tetra_root(n) = iset2
				   tetra_pointer(n) = 
     1					pocket_root(iset2)
				   pocket_root(iset2)=
     1					pocket_root(iset1)
				   pocket_count(iset2)=
     1					pocket_count(iset2)+
     2					pocket_count(iset1)
				   pocket_count(iset1) = 0
c
				endif
c
800			continue
c
900		continue
c
1000	continue
c
	do 1050 i = 1,ntetra
		tetra_root(i) = 0
1050	continue
c
	k = 0
	do 1200 i = 1,npockets
c
		if(pocket_count(i).eq.0) goto 1200
		k = k + 1
		l = pocket_root(i)
		do 1100 j = 1,pocket_count(i)
			tetra_root(l) = k
			l = tetra_pointer(l)
1100		continue
1200	continue
c
	npockets = k
c
	return
	end
c	Hpsort_int.f	Version 1 6/3/1995	Patrice Koehl
c
c	This subroutine rearranges an array of integers in ascending order
c
	subroutine hpsort_int(ra,n)
c
	integer	n,i,ir,j,l
c
	integer	rra
	integer	ra(n)
c
	if(n.lt.2) return
c
	l = n/2 + 1
	ir = n
100	continue
	if(l.gt.1) then
		l = l - 1
		rra = ra(l)
	else
		rra = ra(ir)
		ra(ir) = ra(1)
		ir = ir -1
		if(ir.eq.1) then
			ra(1) = rra
			return
		endif
	endif
	i = l
	j = l+l
200	if(j.le.ir) then
		if(j.lt.ir) then
			if(ra(j).lt.ra(j+1)) j = j + 1
		endif
		if(rra.lt.ra(j)) then
			ra(i) = ra(j)
			i = j
			j = j + j
		else
			j = ir + 1
		endif
		goto 200
	endif
	ra(i) = rra
	goto 100
c
	end
c	Measure_pocket.f	Version 1 8/4/1998	Patrice Koehl
c
c	Based on the alpha shape, this subroutine computes the surface
c	and volume of each pocket (cavity) in the space filling
c	diagram
c
	subroutine measure_pocket(surf_pocket,vol_pocket,pocket_nvertex,
     1	pocket_vertex)
c
	include 'pocket.h'
c
c	Misc. variables
c
	integer	i,j,k,l,m
	integer	ipocket,idx,pair1,pair2,pair3,pair4,pair5,pair6
	integer	trig1,trig2,trig3,trig4
	integer ntetra,ntrig,nedge
	integer npoints,nvertex,option
c
c	Information on voids
c
	integer	npockets,pocket_nvertex
	integer	tetra_pocket(ntetra_max)
	integer	pocket_vertex(npointmax)
	integer	vertex(npointmax)
	integer	index(nvoid_max),rank(nvoid_max)
c
c	Information on tetrahedra
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
	integer tetra(4,ntetra_max)
	integer tetra_neighbour(4,ntetra_max)
	integer tetra_link(4,ntetra_max)
c
c	Information on triangles
c
	integer trig(3,ntrig_max)
	integer trig_link(3,ntrig_max)
	integer	link_trig(2,ntrig_max)
	integer*1 trig_status(ntrig_max)
	integer*1 trig_coef(ntrig_max)
	integer*1 nlink_trig(ntrig_max)
c
c	Information on edges
c
	integer edge(2,nedge_max)
	integer*1 edge_status(nedge_max)
c
c	Information on vertices
c
	integer redinfo(npointmax)
c
	real*8	tetra_volume
	real*8	ra,rb,rc,rd,ra2,rb2,rc2,rd2
	real*8	wa,wb,wc,wd
	real*8	surfa,surfb,surfc
	real*8	vola,volb,volc,vol,coef
	real*8	dist1,dist2,dist3,dist4,dist5,dist6
	real*8	d2_1,d2_2,d2_3,d2_4,d2_5,d2_6
	real*8	ang1,ang2,ang3,ang4,ang5,ang6
	real*8	a(3),b(3),c(3),d(3)
	real*8	surf_pocket(nvoid_max),vol_pocket(nvoid_max)
	real*8	temp(nvoid_max)
	real*8	dsurfa(3,2),dsurfb(3,2)
	real*8  dvola(3,2),dvolb(3,2)
c
	real*8  radius(npointmax),radius2(npointmax)
	real*8	weight(npointmax)
	real*8  coord(3*npointmax)
c
	common	/xyz_vertex/	coord,radius,weight
	common	/tetra_zone/	ntetra,tetra,tetra_neighbour
	common	/tetra_stat/	tetra_status,tetra_orient,
     1				tetra_nindex
	common	/trig_zone/	ntrig,trig
	common	/trig_stat/	trig_status,trig_coef
	common	/edge_zone/	nedge,edge
	common	/edge_stat/	edge_status
	common	/vertex_zone/	npoints,nvertex,redinfo
	common	/links/		tetra_link,trig_link
	common /trig_link/	nlink_trig,link_trig
	common /pockets/	npockets,tetra_pocket
c	
	include 'constants.h'
c
	option = 0
c
	do 100 i = 1,npoints
		radius2(i) = radius(i)*radius(i)
100	continue
c
c
	do 400 ipocket = 1,npockets
c
	    surf_pocket(ipocket) = 0
	    vol_pocket(ipocket) = 0
c
c
	    do 300 idx = 1,ntetra
c
		if(tetra_pocket(idx).ne.ipocket) goto 300
c
		i = tetra(1,idx)
		j = tetra(2,idx)
		k = tetra(3,idx)
		l = tetra(4,idx)
c
		ra = radius(i)
		rb = radius(j)
		rc = radius(k)
		rd = radius(l)
c
		ra2 = radius2(i)
		rb2 = radius2(j)
		rc2 = radius2(k)
		rd2 = radius2(l)
c
		wa = 0.5*weight(i)
		wb = 0.5*weight(j)
		wc = 0.5*weight(k)
		wd = 0.5*weight(l)
c
		trig1 = tetra_link(4,idx)
		trig2 = tetra_link(3,idx)
		trig3 = tetra_link(2,idx)
		trig4 = tetra_link(1,idx)
c
		pair1 = trig_link(3,trig1)
		pair2 = trig_link(2,trig1)
		pair4 = trig_link(1,trig1)
		pair3 = trig_link(2,trig2)
		pair5 = trig_link(1,trig2)
		pair6 = trig_link(1,trig3)
c
		do 200 m = 1,3
			a(m) = coord(3*(i-1)+m)
			b(m) = coord(3*(j-1)+m)
			c(m) = coord(3*(k-1)+m)
			d(m) = coord(3*(l-1)+m)
200		continue
c
		call tetra_6dihed(a,b,c,d,ang1,ang2,ang4,ang3,
     1		ang5,ang6)
c
		vol = tetra_volume(a,b,c,d)
c
		vol_pocket(ipocket) = vol_pocket(ipocket) + vol
c
c	First check all vertices of the tetrahedron
c
		if(redinfo(i).eq.0) then
c
			vertex(i) = ipocket
			coef = 0.5*(ang1+ang2+ang3-0.5)
			surfa = 4*pi*radius2(i)*coef
			vola = surfa*radius(i)/3
c
			surf_pocket(ipocket)=surf_pocket(ipocket)+surfa
			vol_pocket(ipocket)=vol_pocket(ipocket)-vola
c
		endif
c
		if(redinfo(j).eq.0) then
c
			vertex(j) = ipocket
			coef = 0.5*(ang1+ang4+ang5-0.5)
			surfa = 4*pi*radius2(j)*coef
			vola = surfa*radius(j)/3
c
			surf_pocket(ipocket)=surf_pocket(ipocket)+surfa
			vol_pocket(ipocket)=vol_pocket(ipocket)-vola
c
		endif
c
		if(redinfo(k).eq.0) then
c
			vertex(k) = ipocket
			coef = 0.5*(ang2+ang4+ang6-0.5)
			surfa = 4*pi*radius2(k)*coef
			vola = surfa*radius(k)/3
c
			surf_pocket(ipocket)=surf_pocket(ipocket)+surfa
			vol_pocket(ipocket)=vol_pocket(ipocket)-vola
c
		endif
c
		if(redinfo(l).eq.0) then
c
			vertex(l) = ipocket
			coef = 0.5*(ang3+ang5+ang6-0.5)
			surfa = 4*pi*radius2(l)*coef
			vola = surfa*radius(l)/3
c
			surf_pocket(ipocket)=surf_pocket(ipocket)+surfa
			vol_pocket(ipocket)=vol_pocket(ipocket)-vola
c
		endif
c
c	Now check all edges
c
		if(edge_status(pair1).eq.1) then
c
			call distance2(coord,i,j,d2_1,ncortot)
			dist1 = sqrt(d2_1)
c
			call twosphere_vol(a,b,ra,ra2,rb,rb2,dist1,
     1			d2_1,surfa,surfb,vola,volb,
     2			dsurfa,dsurfb,dvola,dvolb,option)
c
			surf_pocket(ipocket)=surf_pocket(ipocket) -
     1			ang1*(surfa+surfb)
			vol_pocket(ipocket) = vol_pocket(ipocket) +
     1			ang1*(vola+volb)
c
		endif
c
		if(edge_status(pair2).eq.1) then
c
			call distance2(coord,i,k,d2_2,ncortot)
			dist2 = sqrt(d2_2)
c
			call twosphere_vol(a,c,ra,ra2,rc,rc2,dist2,
     1			d2_2,surfa,surfb,vola,volb,
     2			dsurfa,dsurfb,dvola,dvolb,option)
c
			surf_pocket(ipocket)=surf_pocket(ipocket) -
     1			ang2*(surfa+surfb)
			vol_pocket(ipocket) = vol_pocket(ipocket) +
     1			ang2*(vola+volb)
c
		endif
c
		if(edge_status(pair3).eq.1) then
c
			call distance2(coord,i,l,d2_3,ncortot)
			dist3 = sqrt(d2_3)
c
			call twosphere_vol(a,d,ra,ra2,rd,rd2,dist3,
     1			d2_3,surfa,surfb,vola,volb,
     2			dsurfa,dsurfb,dvola,dvolb,option)
c
			surf_pocket(ipocket)=surf_pocket(ipocket) -
     1			ang3*(surfa+surfb)
			vol_pocket(ipocket) = vol_pocket(ipocket) +
     1			ang3*(vola+volb)
c
		endif
c
		if(edge_status(pair4).eq.1) then
c
			call distance2(coord,j,k,d2_4,ncortot)
			dist4 = sqrt(d2_4)
c
			call twosphere_vol(b,c,rb,rb2,rc,rc2,dist4,
     1			d2_4,surfa,surfb,vola,volb,
     2			dsurfa,dsurfb,dvola,dvolb,option)
c
			surf_pocket(ipocket)=surf_pocket(ipocket) -
     1			ang4*(surfa+surfb)
			vol_pocket(ipocket) = vol_pocket(ipocket) +
     1			ang4*(vola+volb)
c
		endif
c
		if(edge_status(pair5).eq.1) then
c
			call distance2(coord,j,l,d2_5,ncortot)
			dist5 = sqrt(d2_5)
			call twosphere_vol(b,d,rb,rb2,rd,rd2,dist5,
     1			d2_5,surfa,surfb,vola,volb,
     2			dsurfa,dsurfb,dvola,dvolb,option)
c
			surf_pocket(ipocket)=surf_pocket(ipocket) -
     1			ang5*(surfa+surfb)
			vol_pocket(ipocket) = vol_pocket(ipocket) +
     1			ang5*(vola+volb)
c
		endif
c
		if(edge_status(pair6).eq.1) then
c
			call distance2(coord,k,l,d2_6,ncortot)
			dist6 = sqrt(d2_6)
c
			call twosphere_vol(c,d,rc,rc2,rd,rd2,dist6,
     1			d2_6,surfa,surfb,vola,volb,
     2			dsurfa,dsurfb,dvola,dvolb,option)
c
			surf_pocket(ipocket)=surf_pocket(ipocket) -
     1			ang6*(surfa+surfb)
			vol_pocket(ipocket) = vol_pocket(ipocket) +
     1			ang6*(vola+volb)
c
		endif
c
c	Finally check faces
c
		if(trig_status(trig1).eq.1) then
c
			call threesphere_vol(a,b,c,ra,rb,rc,ra2,
     1			rb2,rc2,wa,wb,wc,dist1,dist2,dist4,
     2			d2_1,d2_2,d2_4,surfa,surfb,surfc,
     3			vola,volb,volc)
c
			surf_pocket(ipocket) = surf_pocket(ipocket)+
     1			0.5*(surfa+surfb+surfc)
			vol_pocket(ipocket) = vol_pocket(ipocket)-
     1			0.5*(vola+volb+volc)
c
		endif
c
		if(trig_status(trig2).eq.1) then
c
			call threesphere_vol(a,b,d,ra,rb,rd,ra2,
     1			rb2,rd2,wa,wb,wd,dist1,dist3,dist5,
     2			d2_1,d2_3,d2_5,surfa,surfb,surfc,
     3			vola,volb,volc)
c
			surf_pocket(ipocket) = surf_pocket(ipocket)+
     1			0.5*(surfa+surfb+surfc)
			vol_pocket(ipocket) = vol_pocket(ipocket)-
     1			0.5*(vola+volb+volc)
c
		endif
c
		if(trig_status(trig3).eq.1) then
c
			call threesphere_vol(a,c,d,ra,rc,rd,ra2,
     1			rc2,rd2,wa,wc,wd,dist2,dist3,dist6,
     2			d2_2,d2_3,d2_6,surfa,surfb,surfc,
     3			vola,volb,volc)
c
			surf_pocket(ipocket) = surf_pocket(ipocket)+
     1			0.5*(surfa+surfb+surfc)
			vol_pocket(ipocket) = vol_pocket(ipocket)-
     1			0.5*(vola+volb+volc)
c
		endif
c
		if(trig_status(trig4).eq.1) then
c
			call threesphere_vol(b,c,d,rb,rc,rd,rb2,
     1			rc2,rd2,wb,wc,wd,dist4,dist5,dist6,
     2			d2_4,d2_5,d2_6,surfa,surfb,surfc,
     3			vola,volb,volc)
c
			surf_pocket(ipocket) = surf_pocket(ipocket)+
     1			0.5*(surfa+surfb+surfc)
			vol_pocket(ipocket) = vol_pocket(ipocket)-
     1			0.5*(vola+volb+volc)
c
		endif
c
300	   continue
c
400	continue
c
	do 500 i = 1,npockets
		if(surf_pocket(i).lt.0) surf_pocket(i) = 0
		if(vol_pocket(i).lt.0) vol_pocket(i) = 0
500	continue
c
	call hpsort_key_down(surf_pocket,index,npockets)
c
	do 600 i = 1,npockets
		temp(i) = vol_pocket(index(i))
600	continue
c
	do 700 i = 1,npockets
		vol_pocket(i) = temp(i)
		rank(index(i)) = i
700	continue
c
	do 800 i = 1,ntetra
c
		j = tetra_pocket(i)
		if(j.ne.0) then
			tetra_pocket(i) = rank(j)
		endif
c
800	continue
c
	pocket_nvertex = 0
c
	do 1200 ipocket = 1,npockets
c
	    pocket_nvertex = pocket_nvertex + 1
	    pocket_vertex(pocket_nvertex) = - ipocket
c
	    do 900 i = 1,npoints
		vertex(i) = 0
900	    continue
c
	    do 1000 idx = 1,ntetra
c
		if(tetra_pocket(idx).ne.ipocket) goto 1000
c
		i = tetra(1,idx)
		j = tetra(2,idx)
		k = tetra(3,idx)
		l = tetra(4,idx)
		if(redinfo(i).eq.0) vertex(i) = ipocket
		if(redinfo(j).eq.0) vertex(j) = ipocket
		if(redinfo(k).eq.0) vertex(k) = ipocket
		if(redinfo(l).eq.0) vertex(l) = ipocket
1000	    continue
c
	    do 1100 i = 1,npoints
		if(vertex(i).ne.0) then
			pocket_nvertex=pocket_nvertex+1
			pocket_vertex(pocket_nvertex)=i
		endif
		vertex(i) = 0
1100	    continue
c
1200	continue
c
	return
	end
c	Hpsort_key_down.f	Version 1 6/3/1995	Patrice Koehl
c
c	This subroutine rearranges an array in descending order, and
c	provide an index of the ranked element
c
	subroutine hpsort_key_down(ra,index,n)
c
	integer	n,i,ir,j,l,idx
	integer	index(n)
c
	real*8	rra
	real*8	ra(n)
c
	do 50 i = 1,n
		index(i) = i
50	continue
c
	if(n.lt.2) return
c
	l = n/2 + 1
	ir = n
100	continue
	if(l.gt.1) then
		l = l - 1
		rra = ra(l)
		idx = l
	else
		rra = ra(ir)
		idx = index(ir)
		ra(ir) = ra(1)
		index(ir) = index(1)
		ir = ir -1
		if(ir.eq.1) then
			ra(1) = rra
			index(1) = idx
			return
		endif
	endif
	i = l
	j = l+l
200	if(j.le.ir) then
		if(j.lt.ir) then
			if(ra(j).gt.ra(j+1)) j = j + 1
		endif
		if(rra.gt.ra(j)) then
			ra(i) = ra(j)
			index(i) = index(j)
			i = j
			j = j + j
		else
			j = ir + 1
		endif
		goto 200
	endif
	ra(i) = rra
	index(i) = idx
	goto 100
c
	end
