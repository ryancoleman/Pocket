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
c	depth.f		Version 1 11/24/2000	Patrice Koehl
c
c	This subroutine defines the depth of each tetrahedron of the
c	Delaunay
c
	subroutine depth
c
	include 'pocket.h'
c
	integer j,k
	integer	ntetra
	integer idx,INF,dmax
	integer	nvisited,ncomplete
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
	integer*1 tetra_visited(ntetra_max)
c
	integer	tetra_depth(ntetra_max)
c
	common  /tetra_zone/    ntetra,tetra,tetra_neighbour
	common  /tetra_stat/    tetra_status,tetra_orient,
     1                          tetra_nindex
	common  /tetra_flux/	tetra_flow
	common  /tet_depth/	tetra_depth
c
c	Loop over all tetrahedra:
c	Define their depth as their index, and then test for "flow"
c	to "infinity"
c
	INF = ntetra + 1
c
	nvisited = 0
c
	do 200 idx = 1,ntetra
c
c		"Dead" tetrahedron are ignored
c
		if(tetra_status(idx).eq.0) then
			tetra_visited(idx) = 1
			nvisited = nvisited + 1
			goto 200
		endif
c
		tetra_depth(idx) = idx
c
		do 100 j = 1,4
c
			k = tetra_neighbour(j,idx)
c
			if(k.eq.0.and.tetra_flow(j,idx).eq.1) then
				tetra_depth(idx) = INF
				tetra_visited(idx) = 1
				nvisited = nvisited + 1
				goto 200
			endif
c
100		continue
c
		tetra_visited(idx) = 0
c
200	continue
c
c	Mimmick a recursive algorithm in C
c
300	continue
c
	if(nvisited.eq.ntetra) goto 600
c
	do 500 idx = ntetra,1,-1
c
		if(tetra_visited(idx).eq.1) goto 500
c
		dmax = tetra_depth(idx)
		ncomplete = 1
c
		do 400 j = 1,4
c
		   k = tetra_neighbour(j,idx)
		   if(tetra_flow(j,idx).eq.1) then
			if(tetra_visited(k).eq.1) then
				if(tetra_depth(k).eq.INF) then
					tetra_depth(idx)=INF
					nvisited = nvisited+1
					tetra_visited(idx) = 1
					goto 500
				endif
				dmax=max(dmax,tetra_depth(k))
			else
				ncomplete = 0
			endif
		   endif
c
400		continue
c
		if(ncomplete.eq.1) then
			tetra_depth(idx) = dmax
			nvisited = nvisited + 1
			tetra_visited(idx) = 1
		endif
c
500	continue
c
	goto 300
c
600	continue
c
	return
	end
