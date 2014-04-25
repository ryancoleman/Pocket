c       alpha_out.f version 1 17 October 2012 Ryan G. Coleman
c
c       this writes out the atoms, edges, triangles and tetrahedra
c       of the alpha complex, useful for other code that can read this in
c       and process the alpha complex data
c
	subroutine alpha_out(label, radius_h2o, radius_sph)
c
	include 'pocket.h'
c
        character label(natot)*54
        real*8  radius_h2o
        real*8  radius_sph(natot)
c
	integer i
	integer	ntetra,ntetraout
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
        integer edge(2,nedge_max)
        integer	trig(3,ntrig_max)
c
        character name*4
        character resname*4,chname*1
        integer ires
c
        common  /vertex_zone/   npoints,nvertex,redinfo
	common  /xyz_vertex/	coord,radius,weight
	common  /tetra_zone/    ntetra,tetra,tetra_neighbour
	common  /tetra_stat/    tetra_status,tetra_orient,
     1                          tetra_nindex
	common  /tetra_flux/	tetra_flow
	common  /gmp_info/	scale,eps
        common  /trig_zone/     ntrig,trig
        common /trig_stat/	trig_status,trig_coef
        common  /edge_zone/     nedge,edge
        common  /edge_stat/     edge_status
        common  /links/         tetra_link,trig_link
        common  /hull/          tetra_hull,trig_hull,edge_hull,
     1                          vertex_hull
c format line for points
1000    format('POINT',x,i10,x,a4,x,a4,x,a1,1x,i4,1x,f8.3,x,f8.3,x,
     1      f8.3,x,f6.3)
c       first write out the points
        do i = 1,npoints
          read(label(i),'(13x,a4,a4,a1,i4)') name,resname,chname,ires
c         we use i-1 as the output index, because this will be read by code
c         that uses 0-indexed arrays not 1-indexed arrays
          write(10,1000) i-1,name,resname,chname,ires,
     1        (coord(3*(i-1)+k),k=1,3),radius_sph(i)-radius_h2o
        enddo
c format line for edges
2000    format('EDGE',x,i10,x,i10,x,i10)
c       now write out the edges
        do i = 1,nedge
c         everything is done - 1 here, since the output is 0-indexed
          write(10,2000) i-1,edge(1,i)-1,edge(2,i)-1     
        enddo
c format line for triangles
3000    format('TRIANGLE',x,i10,x,i10,x,i10,x,i10)
c       write the triangles
        do i = 1,ntrig
c         everything is done - 1 here, since the output is 0-indexed
          write(10,3000) i-1,trig(1,i)-1,trig(2,i)-1,trig(3,i)-1     
        enddo
c format line for tetrahedra
4000    format('TETRA',x,i10,x,i10,x,i10,x,i10,x,i10)
c       write the tetrahedra
        ntetraout = 0
        do i = 1,ntetra
          if (tetra_status(i) .eq. 1) then !this tetrahedra is okay
c           everything is done - 1 here, since the output is 0-indexed
            write(10,4000) ntetraout,tetra(1,i)-1,tetra(2,i)-1,
     1          tetra(3,i)-1,tetra(4,i)-1
            ntetraout = ntetraout + 1
          endif
        enddo
c
	return
	end subroutine alpha_out
