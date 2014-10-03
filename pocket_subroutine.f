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
c	Pocket.f		Version 1 8/4/1998	Patrice Koehl
c
c	This program finds the pockets in a protein, and compute
c	their surface and volume
c
c	Radii of each atom can be arbitrarily defined
c
c	Only polar hydrogens are included
c
        subroutine pocket_sub

	include 'pocket.h'
c
	integer	i,j,k,isurf,nchain,idx,idx1,keep_h
	integer	nseq,natom,nsphere,ntype,nat
	integer	itest,npocket,nmouth
	integer natom_pocket,ipocket,imouth,iskip
	integer	natom_mouth
	integer	nswap,option,lout
	integer	itype(nrestot),listatom(nrestot),chain(nrestot)
	integer	occup(natot)
	integer polartype(nresdef,20),spheretype(nresdef,20)
	integer solvtype(nresdef,20)
	integer	atom_pocket(npointmax)
	integer	atom_mouth(npointmax)
	integer	tetra_pocket(ntetra_max)
	integer	trig_mouth(ntrig_max)
	integer	mouth_pocket(nmouth_max)
	integer	pocket_nmouth(nvoid_max)
	integer	indx(nmouth_max),irank(nmouth_max)
c
	real*8	alpha
	real*8	radius_h2o
	real*8	surftot_pocket,voltot_pocket
	real*8	surftot,voltot,Ssolv,Vsolv
	real*8	coord(ncortot),radprot(natot)
	real*8	coord_sph(ncortot),radius_sph(natot)
	real*8	aspv(natot),coefasp(natot)
	real*8	radii(nresdef,20),aspval(nresdef,20)
	real*8	asp(naspmax)
	real*8	surf_pocket(nvoid_max),vol_pocket(nvoid_max)
	real*8	surf_mouth(nmouth_max)
	real*8	surf(natot),vol(natot)
	real*8  dsurf(3,natot),dvol(3,natot)
c
	character	seq(nrestot)*4,atom(natot)*4
	character	label(natot)*54,label_pdb(natot)*54
	character	chname(10)*1
	character	fname1*50,fname2*50,protname*10,fname3*50
	character	fname4*50,fname5*50,fname6*50,fname7*50
	character	fname8*50,fname9*50,fname10*50
c
        common /protein/ 	nseq,itype,chain,seq
        common /xyz/     	coord,occup,natom,listatom
        common /names/   	atom
	common /polarity/  	polartype,spheretype,solvtype
	common /aspval/     	asp
	common /pockets/	npocket,tetra_pocket
	common /mouths/		nmouth,mouth_pocket,trig_mouth
c
1       format(a)
2       format(1x,'PDB file for all atoms of the protein        : ',$)
3       format(1x,'Generic name for outputs                     : ',$)
5	format(1x,'Enter radius of water probe                  : ',$)
8	format(1x,'Parameter file for VdW radii                 : ',$)
10	format(1x,'Name of generic input file                   : ',$)
16	format(/,'REMARK',28x,'X',7x,'Y',7x,'Z',7x,'Rvdw',3x,'Probe',
     &  8x,'Surf',9x,'Vol',2x,' dSurft/dX',2x,' dSurft/dY',2x,
     &  ' dSurft/dZ',2x,'  dVolt/dX',2x,'  dVolt/dY',2x,'  dVolt/dZ')
22	format(' Id ',2x,'N_mth',2x,'     Surface',
     &  2x,'      Volume')
23	format(i4,2x,i5,2x,2(f12.4,2x),2x,i4)
24	format(/,'Tot ',9x,2(f12.4,2x),2x,i4)
26	format(a54,14x,i4,2x,a4)
28	format('@kinemage 1',/)
29	format('@caption',/)
30	format('@onewidth',/,'@whitebkg',/,/)
31	format(a54,2x,f6.2,2x,f6.2,2x,8(f10.4,2x))
32	format(' Id ',2x,'Pocket',2x,'    Surface')
33	format(i4,3x,i4,2x,f12.4)
34	format(a54,14x,i4,2x,a4)
c
        write(6,*) ' '
c
	ntype = 1
c
	if(ntype.eq.1) then
c
        	write(6,2)
        	read(5,1) fname1
c
        	write(6,8)
        	read(5,1) fname3
c
		write(6,3)
		read(5,1) protname
c
		isurf = 1
c
		if(isurf.eq.1) then
			write(6,5)
			read(5,*) radius_h2o
		else
			radius_h2o = 0.0
		endif
c
		keep_h = 0
c
	else
		write(6,10)
		read(5,1) fname1
	endif
c
c      	write(6,3)
c       read(5,1) fname2
c
c	write(6,13)
c	read(5,1) fname4
c
c	If file is a PDB file, define topology, read atomic radii,
c	coordinates of all atoms, and define spheres with their radii
c
	if(ntype.eq.1) then
c
c		Read amino acid topologies and parameters
c
		call topology_new
		call readrad(fname3,radii,aspval,itest)
c
c       	Read pdb file to get sequence and coordinates 
c		information
c
		call defchain(fname1,nchain,chname)
        	call extractseq(fname1)
        	call prepare
        	call readpdb(fname1,keep_h,label_pdb)
c
c		Define radius of each atom
c
		nat = 0
		do 200 i = 1, nseq
			do 100 j = 1,listatom(i)
				nat = nat + 1
				radprot(nat)=radii(itype(i),j)+
     1					radius_h2o
				if(itest.eq.0) then
					aspv(nat) = asp(
     1					solvtype(itype(i),j))
				else
					aspv(nat) = 
     1					aspval(itype(i),j)
				endif
100			continue
200		continue
c
c		Define coordinates of the center of the spheres 
c		considered, as well as VdW radii
c
		nsphere = 0
		do 400 i = 1,natom
c
			if(keep_h.eq.0.and.atom(i)(1:1).eq.'H') 
     1			goto 400
			if(occup(i).eq.0) goto 400
c
			nsphere = nsphere + 1
			do 300 j = 1,3
				coord_sph(3*(nsphere-1)+j)=
     1					coord(3*(i-1)+j)
300			continue
c
			radius_sph(nsphere) = radprot(i)
			coefasp(nsphere)    = aspv(i)
			coefasp(nsphere)    = 1.
			label(nsphere)	    = label_pdb(i)
c
400		continue
c
	else
c
c	The file contains a generic collection of points, with the 
c	number of points on the first line, and the coordinates + 
c	radius of each point on the following lines
c
		open(unit=1,file=fname1,status='old')
		read(1,*) nsphere
		do 500 i = 1,nsphere
			read(1,*) (coord_sph(3*(i-1)+j),j=1,3),
     1				radius_sph(i)
			coefasp(i) = i
c			coefasp(i) = 1.
500		continue
		close(unit=1)
c
	endif
c
c	Now calculate volume and surface
c
c	1. Perform the regular triangulation
c
	call regular3D(coord_sph,radius_sph,nsphere)
c
c	2. Define all simplices of the regular triangulation
c
	call define_triangles
	call define_edges
c
c	Define flow
c
	call delaunay_flow
	call depth
c
c	4. Compute the Alpha complex for fixed value of alpha
c
	alpha = 0.
	call alpha_fixed(alpha)
c
c	5. Compute surface and volume of the protein
c
	option = 1
	call measure_vol(coefasp,Ssolv,Vsolv,surftot,surf,voltot,vol,
     1		dsurf,dvol,option)
c
c	6. Now compute pockets
c
	call find_pockets
	call measure_pocket(surf_pocket,vol_pocket,natom_pocket,
     1	atom_pocket)
c
c	7. Now look at mouths of pockets
c
	call find_mouth
	call measure_mouth(surf_mouth,natom_mouth,atom_mouth)
c
	do 600 i = 1,npocket
		pocket_nmouth(i) = 0
600	continue
c
	do 700 i = 1,nmouth
		pocket_nmouth(mouth_pocket(i))=
     1			pocket_nmouth(mouth_pocket(i))+1
700	continue
c
c	Now store the results
c
c	Get protein name
c
	idx1 = index(protname,' ') -1
	fname4=protname(1:idx1)//".pocket"
	fname2=protname(1:idx1)//".info"
	fname5=protname(1:idx1)//".rasmol"
	fname6=protname(1:idx1)//".kin"
c	fname7=protname(1:idx1)//".vol"
	fname8=protname(1:idx1)//".atvol"
	fname9=protname(1:idx1)//".mouth"
	fname10=protname(1:idx1)//".alpha"

c
	open(unit=2,file=fname2,status='unknown')
	open(unit=4,file=fname4,status='unknown')
c	open(unit=7,file=fname7,status='unknown')
	open(unit=8,file=fname8,status='unknown')
	open(unit=9,file=fname9,status='unknown')
	open(unit=10,file=fname10,status='unknown')
c
	lout = 11
	open(unit=lout,file=fname6,status='unknown')
c
c	First write surf+volume:
c	in summary file (.vol), and per atom (.atvol)
c
c	write(7,27)
c	write(7,7) protname(1:idx1),nseq-2,nsphere,radius_h2o,
c     1		surftot,voltot
c	close(unit=7)
c
	write(8,16)
	do 800 i = 1,nsphere
		write(8,31) label(i),radius_sph(i)-radius_h2o,
     1		radius_h2o,surf(i),vol(i),(dsurf(j,i),j=1,3),
     2		(dvol(j,i),j=1,3)
800	continue
	close(unit=8)
c
	write(6,*) ' '
	write(6,*) 'Protein             : ',protname
	write(6,*) 'Probe radius        : ',sngl(radius_h2o)
	write(6,*) 'Number of atoms     : ',nsphere
	write(6,*) 'Total Surface Area  : ',sngl(surftot)
	write(6,*) 'Total Volume        : ',sngl(voltot)
	write(6,*) ' '
c
	write(2,*) ' '
	write(2,*) 'Protein             : ',protname
	write(2,*) 'Probe radius        : ',sngl(radius_h2o)
	write(2,*) 'Number of atoms     : ',nsphere
	write(2,*) 'Total Surface Area  : ',sngl(surftot)
	write(2,*) 'Total Volume        : ',sngl(voltot)
	write(2,*) ' '
c
	write(lout,28)
	write(lout,29)
	write(lout,*) 'Protein 		   : ',protname
	write(lout,*) 'Probe radius	   : ',sngl(radius_h2o)
	write(lout,*) 'Number of atoms     : ',nsphere
	write(lout,*) 'Total Surface Area  : ',sngl(surftot)
	write(lout,*) 'Total Volume	   : ',sngl(voltot)
	write(lout,*) ' '
c
c	Now write pocket information
c
	surftot_pocket = 0
	voltot_pocket = 0
c
	write(6,*) ' '
	write(6,*) 'Pockets : '
	write(2,*) ' '
	write(2,*) 'Pockets : '
	write(lout,*) 'Pockets : '
	write(6,*) ' '
	write(2,*) ' '
	write(lout,*) ' '
	write(2,22)
	write(6,22)
	write(lout,22)
	do 900 i = 1,npocket
		write(2,23) i,pocket_nmouth(i),surf_pocket(i),
     1		vol_pocket(i)
		write(lout,23) i,pocket_nmouth(i),surf_pocket(i),
     1		vol_pocket(i)
		write(6,23) i,pocket_nmouth(i),surf_pocket(i),
     1		vol_pocket(i)
		surftot_pocket = surftot_pocket + surf_pocket(i)
		voltot_pocket = voltot_pocket + vol_pocket(i)
900	continue
	write(2,24) surftot_pocket,voltot_pocket
	write(6,24) surftot_pocket,voltot_pocket
c	write(lout,24) surftot_pocket,voltot_pocket
c
	write(6,*) ' '
	write(6,*) 'Mouths : '
	write(2,*) ' '
	write(2,*) 'Mouths : '
	write(lout,*) ' '
	write(lout,*) 'Mouths : '
	write(6,*) ' '
	write(2,*) ' '
	write(lout,*) ' '
	write(2,32)
	write(6,32)
	write(lout,32)
c
	call isort_indx(mouth_pocket,indx,nswap,nmouth)
c
	k = 0
	do 1100 i = 1,nmouth
		j = indx(i)
		if(surf_mouth(j).ne.0) then
			k = k+1
			write(2,33) k,mouth_pocket(i),
     1				surf_mouth(j)
			write(lout,33) k,mouth_pocket(i),
     1				surf_mouth(j)
			write(6,33) k,mouth_pocket(i),
     1				surf_mouth(j)
			irank(i) = k
		else
			irank(i) = 0
		endif
1100	continue
c
	ipocket = 0
	do 1200 i = 1,natom_pocket
		if(atom_pocket(i).lt.0) then
			ipocket = -atom_pocket(i)
			goto 1200
		endif
		j = atom_pocket(i)
		write(4,26) label(j),ipocket,'POC '
1200	continue
c
	imouth = 0
	do 1300 i = 1,natom_mouth
		if(atom_mouth(i).lt.0) then
			imouth = -atom_mouth(i)
			if(irank(imouth).eq.0) then
				iskip = 1
			else
				iskip = 0
			endif
			goto 1300
		endif
		if(iskip.eq.0) then
			j = atom_mouth(i)
			write(9,34) label(j),irank(imouth),
     1				'MTH '
		endif
1300	continue
c
	close(unit=2)
	close(unit=4)
	close(unit=9)

!writes out a bunch of stuff to unit=10, the a-complex
        call alpha_out(label, radius_h2o, radius_sph) 

	close(unit=10)
c
	call rasmol_script(fname5,label,natom_pocket,atom_pocket)
c
	write(lout,*) ' '
	write(lout,30)
c
	call kin_script(lout,label,label_pdb,radii,radius_h2o,
     1			natom_pocket,atom_pocket)
c
	close(unit=lout)
c
	end subroutine
