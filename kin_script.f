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
c	kin_script.f		Version 1 5/14/2002	Patrice Koehl
c
c	This program generates a KineMage script to view a molecule,
c	visualize the Delaunay triangulation based on its atoms,
c	as well as its pockets.
c
	subroutine kin_script(lout,label,label_pdb,radii,radius_h2o,
     1			natom_pocket,atom_pocket)
c
	include 'pocket.h'
c
	integer	lout
	integer	natom_pocket
	integer	atom_pocket(natot)
c
	character label(natot)*54
	character label_pdb(natot)*54
	character tag(natot)*14
c
	real*8	radius_h2o
	real*8	radii(nresdef,20)
c
c	Start protein
c
	call kin_protein(lout,label_pdb)
c
c	Add CPK model (solvated, and non solvated)
c
	call kin_cpk(lout,label_pdb,radii,radius_h2o)
c
c	Now build Delaunay
c
	call kin_delaunay(lout,label,tag)
c
c	Now add Pocket in CPK
c
	call kin_pocket(lout,tag,radius_h2o,
     1			natom_pocket,atom_pocket)
c
	return
	end
c
c	kin_protein.f		Version 1 5/14/2002	Patrice Koehl
c
c	This program generates a kin script to view a molecule in
c	PDB format
c
	subroutine kin_protein(lout,label)
c
	include 'pocket.h'
c
	real*8	coord(ncortot)
c
	integer	nat,i,j,k,natom,nseq,ires
	integer	lout
	integer	itype(nrestot),listatom(nrestot)
	integer	chain(nrestot),occup(natot)
c
	character atom(natot)*4,seq(nrestot)*4
	character name*4
	character resname*4,chname*1
	character label(natot)*54
	character name_lc*4,res_lc*4,ch_lc*4
	character save*96
	character to_lower4*4,to_lower1
c
	common /protein/ nseq,itype,chain,seq
	common /xyz/     coord,occup,natom,listatom
	common /names/   atom
c
1	format('@group {Protein} ',/)
2	format('@subgroup {    Mainchain} dominant')
3	format('@subgroup {    Sidechain} dominant')
5	format('@vectorlist {Main} color=gray')
6	format(21x,a1,i4)
7	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
8	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
9	format(a)
10	format('@vectorlist {Side} color=gray')
c
c	Start protein
c
	write(lout,1)
	write(lout,2)
	write(lout,5)
c
c	Write Mainchain
c
	nat = listatom(1)
c
	do 200 i = 2,nseq-1
		resname = seq(i)
		if(i.eq.1) resname = seq(i+1)
		if(i.eq.nseq) resname = seq(i-1)
		ires = i
		if(i.gt.1) ires = ires-1
		if(i.eq.nseq) ires = ires -1
		res_lc = to_lower4(resname)
		do 100 j = 1,listatom(i)
			nat = nat + 1
			if(j.gt.6) goto 100
			if(j.eq.2.or.j.eq.4) goto 100
			if(occup(nat).eq.0) goto 100
			read(label(nat),6) chname,ires
			name = atom(nat)
			name_lc = to_lower4(name)
			ch_lc = to_lower1(chname)
			if(i.eq.2.and.j.eq.1) then
				write(lout,7) name_lc, res_lc,ch_lc,
     1				ires,'P',(coord(3*nat-3+k),k=1,3),
     2				name_lc, res_lc,ch_lc,ires,'L',
     3				(coord(3*nat-3+k),k=1,3)
			else
				write(lout,8) name_lc, res_lc,ch_lc,
     1				ires,'L',(coord(3*nat-3+k),k=1,3)
				if(j.eq.5) then
				write(save,7) name_lc, res_lc,ch_lc,
     1				ires,'P',(coord(3*nat-3+k),k=1,3),
     2				name_lc, res_lc,ch_lc,ires,'L',
     3				(coord(3*nat-3+k),k=1,3)
				endif
			endif
			if(j.eq.6.and.occup(nat-1).eq.1) 
     1				write(lout,9) save
100		continue
200	continue
c
c	Now write sidechains
c
	write(lout,*) ' '
	write(lout,3)
	write(lout,10)
c
	nat = listatom(1)
c
	do 400 i = 2,nseq-1
		if(itype(i).eq.1) then
			nat = nat + listatom(i)
			goto 400
		endif
		resname = seq(i)
		if(i.eq.1) resname = seq(i+1)
		if(i.eq.nseq) resname = seq(i-1)
		ires = i
		if(i.gt.1) ires = ires-1
		if(i.eq.nseq) ires = ires -1
		res_lc = to_lower4(resname)
c
c	Goto CA
c
		nat = nat + 3
		read(label(nat),6) chname,ires
		name = atom(nat)
		name_lc = to_lower4(name)
		ch_lc = to_lower1(chname)
		write(lout,7) name_lc, res_lc,ch_lc,
     1		ires,'P',(coord(3*nat-3+k),k=1,3),
     2		name_lc, res_lc,ch_lc,ires,'L',
     3		(coord(3*nat-3+k),k=1,3)
c
c	Connect CB
c
		nat = nat + 1
		name = atom(nat)
		name_lc = to_lower4(name)
		write(lout,8) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat-3+k),k=1,3)
c
c	Now add rest of the sidechain
c
		nat = nat + 2
c
		call kin_side(lout,itype(i),coord,atom,nat,
     1		res_lc,ch_lc,ires,ncortot)
c
		nat = nat + listatom(i) - 6
c
400	continue
c
	close(unit=1)
c
	return
	end
c
c	to_lower4.f	Version 1 5/14/2002	Patrice Koehl
c
c	This function converts a string of 4 from upper-case to
c	lower-case
c
	function to_lower4(string)
c
	integer	i,j
c
	character*4	string,to_lower4
c
	do 100 i = 1,4
c
		j = ichar(string(i:i))
		if(j.ge.65.and.j.le.90) j = j + 32
		to_lower4(i:i)=char(j)
c
100	continue
c
	return
	end
c
c	to_lower1.f	Version 1 5/14/2002	Patrice Koehl
c
c	This function converts a string of 1 from upper-case to
c	lower-case
c
	function to_lower1(string)
c
	integer	j
c
	character*1	string,to_lower1
c
	j = ichar(string)
	if(j.ge.65.and.j.le.90) then
		j = j+32
	endif
	to_lower1=char(j)
c
	return
	end
c
c	kin_side.f	Version 1 5/14/2002	Patrice Koehl
c
c	This routine generates the bonds for sidechains in proteins
c
	subroutine kin_side(lout,itype,coord,atom,nat,
     1		res_lc,ch_lc,ires,natot)
c
	integer	lout,ires,natot,itype,nat
c
	real*8	coord(3*natot)
c
	character	atom(natot)*4
	character	res_lc*4,ch_lc*1
c
	if(itype.eq.1.or.itype.eq.2) then
		return
	elseif(itype.eq.3) then
		call val_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
	elseif(itype.eq.4) then
		call ile_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
	elseif(itype.eq.5) then
		call leu_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
	elseif(itype.eq.6) then
		call phe_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
	elseif(itype.eq.7) then
		call pro_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
	elseif(itype.eq.8) then
		call met_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
	elseif(itype.eq.9) then
		call trp_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
	elseif(itype.eq.10) then
		call cys_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
	elseif(itype.eq.11) then
		call ser_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
	elseif(itype.eq.12) then
		call thr_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
	elseif(itype.eq.13) then
		call asn_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
	elseif(itype.eq.14) then
		call gln_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
	elseif(itype.eq.15) then
		call tyr_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
	elseif(itype.eq.16) then
		call his_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
	elseif(itype.eq.17) then
		call asp_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
	elseif(itype.eq.18) then
		call glu_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
	elseif(itype.eq.19) then
		call lys_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
	elseif(itype.eq.20) then
		call arg_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
	endif
c
	return
	end
c
c	Val_kin.f	Version 1 5/14/2002	Patrice Koehl
c

	subroutine val_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
c
	integer nat,ires,natot,lout,k
	integer	nat_local
c
	real*8	coord(3*natot)
c
	character atom(natot)*4
	character res_lc*4,ch_lc*1
	character name*4,name_lc*4
	character to_lower4*4
c
1	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
2	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
c
c	nat is at the level of the atom O of the residue
c
	nat_local = nat+1
c
c	connect CB-CG1
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	Move back to CB
c
	nat_local = nat - 2
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,1) name_lc, res_lc,ch_lc,
     1		ires,'P',(coord(3*nat_local-3+k),k=1,3),
     2		name_lc, res_lc,ch_lc,
     3		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	Now connect CG2
c
	nat_local = nat + 2
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
	return
	end
c
c	ile_kin.f	Version 1 5/14/2002	Patrice Koehl
c

	subroutine ile_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
c
	integer nat,ires,natot,lout,k
	integer	nat_local
c
	real*8	coord(3*natot)
c
	character atom(natot)*4
	character res_lc*4,ch_lc*1
	character name*4,name_lc*4
	character to_lower4*4
c
1	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
2	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
c
c	nat is at the level of the atom O of the residue
c
	nat_local = nat+1
c
c	connect CB-CG1
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CG1-CD1
c
	nat_local = nat+3
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	Move back to CB
c
	nat_local = nat - 2
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,1) name_lc, res_lc,ch_lc,
     1		ires,'P',(coord(3*nat_local-3+k),k=1,3),
     2		name_lc, res_lc,ch_lc,
     3		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	Now connect CG2
c
	nat_local = nat + 2
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
	return
	end
c
c	leu_kin.f	Version 1 5/14/2002	Patrice Koehl
c

	subroutine leu_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
c
	integer nat,ires,natot,lout,k
	integer	nat_local
c
	real*8	coord(3*natot)
c
	character atom(natot)*4
	character res_lc*4,ch_lc*1
	character name*4,name_lc*4
	character to_lower4*4
c
1	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
2	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
c
c	nat is at the level of the atom O of the residue
c
	nat_local = nat+1
c
c	connect CB-CG1
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CG1-CD1
c
	nat_local = nat+2
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	Move back to CG1
c
	nat_local = nat + 1
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,1) name_lc, res_lc,ch_lc,
     1		ires,'P',(coord(3*nat_local-3+k),k=1,3),
     2		name_lc, res_lc,ch_lc,
     3		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	Now connect CD2
c
	nat_local = nat + 3
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
	return
	end
c
c	phe_kin.f	Version 1 5/14/2002	Patrice Koehl
c

	subroutine phe_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
c
	integer nat,ires,natot,lout,k
	integer	nat_local
c
	real*8	coord(3*natot)
c
	character atom(natot)*4
	character res_lc*4,ch_lc*1
	character name*4,name_lc*4
	character to_lower4*4
c
1	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
2	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
c
c	nat is at the level of the atom O of the residue
c
	nat_local = nat+1
c
c	connect CB-CG
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CG1-CD1
c
	nat_local = nat+2
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CD1-CE1
c
	nat_local = nat+4
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CE1-CZ
c
	nat_local = nat+6
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CZ-CE2
c
	nat_local = nat+5
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CE2-CD2
c
	nat_local = nat+3
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CD2-CG
c
	nat_local = nat+1
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
	return
	end
c
c	pro_kin.f	Version 1 5/14/2002	Patrice Koehl
c
	subroutine pro_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
c
	integer nat,ires,natot,lout,k
	integer	nat_local
c
	real*8	coord(3*natot)
c
	character atom(natot)*4
	character res_lc*4,ch_lc*1
	character name*4,name_lc*4
	character to_lower4*4
c
1	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
2	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
c
c	nat is at the level of the atom O of the residue
c
c
c	connect CB-CG
c
	nat_local = nat+1
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CG1-CD1
c
	nat_local = nat+2
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CD1-N
c
	nat_local = nat-5
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
	return
	end
c
c	met_kin.f	Version 1 5/14/2002	Patrice Koehl
c
	subroutine met_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
c
	integer nat,ires,natot,lout,k
	integer	nat_local
c
	real*8	coord(3*natot)
c
	character atom(natot)*4
	character res_lc*4,ch_lc*1
	character name*4,name_lc*4
	character to_lower4*4
c
1	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
2	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
c
c	nat is at the level of the atom O of the residue
c
c
c	connect CB-CG
c
	nat_local = nat+1
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CG-SD
c
	nat_local = nat+2
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect SD-CE
c
	nat_local = nat+3
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
	return
	end
c
c	trp_kin.f	Version 1 5/14/2002	Patrice Koehl
c
	subroutine trp_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
c
	integer nat,ires,natot,lout,k
	integer	nat_local
c
	real*8	coord(3*natot)
c
	character atom(natot)*4
	character res_lc*4,ch_lc*1
	character name*4,name_lc*4
	character to_lower4*4
c
1	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
2	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
c
c	nat is at the level of the atom O of the residue
c
	nat_local = nat+1
c
c	connect CB-CG
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CG-CD2
c
	nat_local = nat+3
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CD2-CE3
c
	nat_local = nat+6
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CE3-CZ3
c
	nat_local = nat+8
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CZ3-CH2
c
	nat_local = nat+9
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CH2-CZ2
c
	nat_local = nat+7
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CZ2-CE2
c
	nat_local = nat+5
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CE2-NE1
c
	nat_local = nat+4
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect NE1-CD1
c
	nat_local = nat+2
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CD1-CG
c
	nat_local = nat+1
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	Move back to CD2
c
	nat_local = nat + 3
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,1) name_lc, res_lc,ch_lc,
     1		ires,'P',(coord(3*nat_local-3+k),k=1,3),
     2		name_lc, res_lc,ch_lc,
     3		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CD2-CE2
c
	nat_local = nat+5
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
	return
	end
c
	subroutine cys_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
c
	integer nat,ires,natot,lout,k
	integer	nat_local
c
	real*8	coord(3*natot)
c
	character atom(natot)*4
	character res_lc*4,ch_lc*1
	character name*4,name_lc*4
	character to_lower4*4
c
1	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
2	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
c
c	nat is at the level of the atom O of the residue
c
c	connect CB-SG
c
	nat_local = nat+1
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
	return
	end
c
	subroutine ser_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
c
	integer nat,ires,natot,lout,k
	integer	nat_local
c
	real*8	coord(3*natot)
c
	character atom(natot)*4
	character res_lc*4,ch_lc*1
	character name*4,name_lc*4
	character to_lower4*4
c
1	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
2	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
c
c	nat is at the level of the atom O of the residue
c
c	connect CB-OG
c
	nat_local = nat+1
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
	return
	end
c
c	thr_kin.f	Version 1 5/14/2002	Patrice Koehl
c

	subroutine thr_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
c
	integer nat,ires,natot,lout,k
	integer	nat_local
c
	real*8	coord(3*natot)
c
	character atom(natot)*4
	character res_lc*4,ch_lc*1
	character name*4,name_lc*4
	character to_lower4*4
c
1	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
2	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
c
c	nat is at the level of the atom O of the residue
c
	nat_local = nat+1
c
c	connect CB-CG2
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	Move back to CB
c
	nat_local = nat - 2
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,1) name_lc, res_lc,ch_lc,
     1		ires,'P',(coord(3*nat_local-3+k),k=1,3),
     2		name_lc, res_lc,ch_lc,
     3		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	Now connect OG1
c
	nat_local = nat + 2
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
	return
	end
c
c	asn_kin.f	Version 1 5/14/2002	Patrice Koehl
c
	subroutine asn_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
c
	integer nat,ires,natot,lout,k
	integer	nat_local
c
	real*8	coord(3*natot)
c
	character atom(natot)*4
	character res_lc*4,ch_lc*1
	character name*4,name_lc*4
	character to_lower4*4
c
1	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
2	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
c
c	nat is at the level of the atom O of the residue
c
	nat_local = nat+1
c
c	connect CB-CG
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CG-OD1
c
	nat_local = nat+2
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	Move back to CG
c
	nat_local = nat + 1
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,1) name_lc, res_lc,ch_lc,
     1		ires,'P',(coord(3*nat_local-3+k),k=1,3),
     2		name_lc, res_lc,ch_lc,
     3		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	Now connect ND2
c
	nat_local = nat + 3
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
	return
	end
c
c	gln_kin.f	Version 1 5/14/2002	Patrice Koehl
c
	subroutine gln_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
c
	integer nat,ires,natot,lout,k
	integer	nat_local
c
	real*8	coord(3*natot)
c
	character atom(natot)*4
	character res_lc*4,ch_lc*1
	character name*4,name_lc*4
	character to_lower4*4
c
1	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
2	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
c
c	nat is at the level of the atom O of the residue
c
	nat_local = nat+1
c
c	connect CB-CG
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CG-CD
c
	nat_local = nat+2
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CD-OE1
c
	nat_local = nat+3
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	Move back to CD
c
	nat_local = nat + 2
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,1) name_lc, res_lc,ch_lc,
     1		ires,'P',(coord(3*nat_local-3+k),k=1,3),
     2		name_lc, res_lc,ch_lc,
     3		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	Now connect NE2
c
	nat_local = nat + 4
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
	return
	end
c
c	tyr_kin.f	Version 1 5/14/2002	Patrice Koehl
c

	subroutine tyr_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
c
	integer nat,ires,natot,lout,k
	integer	nat_local
c
	real*8	coord(3*natot)
c
	character atom(natot)*4
	character res_lc*4,ch_lc*1
	character name*4,name_lc*4
	character to_lower4*4
c
1	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
2	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
c
c	nat is at the level of the atom O of the residue
c
	nat_local = nat+1
c
c	connect CB-CG
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CG1-CD1
c
	nat_local = nat+2
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CD1-CE1
c
	nat_local = nat+4
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CE1-CZ
c
	nat_local = nat+6
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CZ-CE2
c
	nat_local = nat+5
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CE2-CD2
c
	nat_local = nat+3
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CD2-CG
c
	nat_local = nat+1
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	Move back to CZ
c
	nat_local = nat+6
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,1) name_lc, res_lc,ch_lc,
     1		ires,'P',(coord(3*nat_local-3+k),k=1,3),
     2		name_lc, res_lc,ch_lc,
     3		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	Connect CZ to OH
c
	nat_local = nat+7
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
	return
	end
c
c	his_kin.f	Version 1 5/14/2002	Patrice Koehl
c

	subroutine his_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
c
	integer nat,ires,natot,lout,k
	integer	nat_local
c
	real*8	coord(3*natot)
c
	character atom(natot)*4
	character res_lc*4,ch_lc*1
	character name*4,name_lc*4
	character to_lower4*4
c
1	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
2	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
c
c	nat is at the level of the atom O of the residue
c
	nat_local = nat+1
c
c	connect CB-CG
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CG1-NE1
c
	nat_local = nat+2
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect NE1-CE1
c
	nat_local = nat+4
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CE1-NE2
c
	nat_local = nat+5
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect NE2-CD2
c
	nat_local = nat+3
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CD2-CG
c
	nat_local = nat+1
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
	return
	end
c
c	asp_kin.f	Version 1 5/14/2002	Patrice Koehl
c
	subroutine asp_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
c
	integer nat,ires,natot,lout,k
	integer	nat_local
c
	real*8	coord(3*natot)
c
	character atom(natot)*4
	character res_lc*4,ch_lc*1
	character name*4,name_lc*4
	character to_lower4*4
c
1	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
2	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
c
c	nat is at the level of the atom O of the residue
c
	nat_local = nat+1
c
c	connect CB-CG
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CG-OD1
c
	nat_local = nat+2
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	Move back to CG
c
	nat_local = nat + 1
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,1) name_lc, res_lc,ch_lc,
     1		ires,'P',(coord(3*nat_local-3+k),k=1,3),
     2		name_lc, res_lc,ch_lc,
     3		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	Now connect OD2
c
	nat_local = nat + 3
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
	return
	end
c
c	glu_kin.f	Version 1 5/14/2002	Patrice Koehl
c
	subroutine glu_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
c
	integer nat,ires,natot,lout,k
	integer	nat_local
c
	real*8	coord(3*natot)
c
	character atom(natot)*4
	character res_lc*4,ch_lc*1
	character name*4,name_lc*4
	character to_lower4*4
c
1	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
2	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
c
c	nat is at the level of the atom O of the residue
c
	nat_local = nat+1
c
c	connect CB-CG
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CG-CD
c
	nat_local = nat+2
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CD-OE1
c
	nat_local = nat+3
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	Move back to CD
c
	nat_local = nat + 2
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,1) name_lc, res_lc,ch_lc,
     1		ires,'P',(coord(3*nat_local-3+k),k=1,3),
     2		name_lc, res_lc,ch_lc,
     3		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	Now connect OE2
c
	nat_local = nat + 4
c
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
	return
	end
c
c	lys_kin.f	Version 1 5/14/2002	Patrice Koehl
c
	subroutine lys_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
c
	integer nat,ires,natot,lout,k
	integer	nat_local
c
	real*8	coord(3*natot)
c
	character atom(natot)*4
	character res_lc*4,ch_lc*1
	character name*4,name_lc*4
	character to_lower4*4
c
1	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
2	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
c
c	nat is at the level of the atom O of the residue
c
c
c	connect CB-CG
c
	nat_local = nat+1
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CG-CD
c
	nat_local = nat+2
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CD-CE
c
	nat_local = nat+3
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CE-NZ
c
	nat_local = nat+4
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
	return
	end
c
c	arg_kin.f	Version 1 5/14/2002	Patrice Koehl
c
	subroutine arg_kin(lout,coord,atom,nat,res_lc,ch_lc,
     1		ires,natot)
c
	integer nat,ires,natot,lout,k
	integer	nat_local
c
	real*8	coord(3*natot)
c
	character atom(natot)*4
	character res_lc*4,ch_lc*1
	character name*4,name_lc*4
	character to_lower4*4
c
1	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
2	format('{ ',a4,a4,a1,1x,i4,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
c
c	nat is at the level of the atom O of the residue
c
c
c	connect CB-CG
c
	nat_local = nat+1
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CG-CD
c
	nat_local = nat+2
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CD-NE
c
	nat_local = nat+3
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect NE-CZ
c
	nat_local = nat+4
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CZ-NH1
c
	nat_local = nat+6
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	return to CZ
c
	nat_local = nat+4
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,1) name_lc, res_lc,ch_lc,
     1		ires,'P',(coord(3*nat_local-3+k),k=1,3),
     2		name_lc, res_lc,ch_lc,
     3		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
c	connect CZ-NH2
c
	nat_local = nat+7
	name = atom(nat_local)
	name_lc = to_lower4(name)
	write(lout,2) name_lc, res_lc,ch_lc,
     1		ires,'L',(coord(3*nat_local-3+k),k=1,3)
c
	return
	end
c
c	kin_pocket.f		Version 1 5/14/2002	Patrice Koehl
c
c	This program generates a KineMage script to view the pockets of a 
c	molecule
c
	subroutine kin_pocket(lout,tag,radius_h2o,natom_pocket,
     1			atom_pocket)
c
	include 'pocket.h'
c
	integer	i,j,k
	integer	lout
	integer	ipocket
	integer	natom_pocket
	integer	atom_pocket(natot)
c
	character tag(natot)*14
c
	real*8	radius_h2o
	real*8	radius(npointmax),weight(npointmax)
	real*8	coord(3*npointmax)
c
	common /xyz_vertex/     coord,radius,weight
c
1	format('@group {Pockets_CPK} ',/)
2	format('@balllist {Ball_Pck} color=green radius = ',f6.3)
4	format('{ ',a14,' }',1x,f8.3,',',f8.3,',',
     1	f8.3)
6	format('@subgroup {    Pocket',i1,'} dominant off')
7	format('@subgroup {    Pocket',i2,'} dominant off')
8	format('@subgroup {    Pocket',i3,'} dominant off')
c
c	Start CPK model
c
	write(lout,1)
c
	ipocket = 0
	do 100 i = 1,natom_pocket
c
		if(atom_pocket(i).lt.0) then
			ipocket = -atom_pocket(i)
			if(ipocket.gt.10) goto 200
			if(ipocket.lt.10) then
				write(lout,6) ipocket
			elseif(ipocket.lt.100) then
				write(lout,7) ipocket
			elseif(ipocket.lt.1000) then
				write(lout,8) ipocket
			endif
			goto 100
		endif
c
		j = atom_pocket(i)
		write(lout,2) radius(j)-radius_h2o
		write(lout,4) tag(j),(coord(3*j-3+k),k=1,3)
c
100	continue
c
200	continue
c
	return
	end
c
c	kin_delaunay.f		Version 1 5/14/2002	Patrice Koehl
c
c	This program generates a KineMage script to view a molecule in
c	PDB format
c
	subroutine kin_delaunay(lout,label,tag)
c
	include 'pocket.h'
c
	real*8	radius(npointmax),weight(npointmax)
	real*8	coord(3*npointmax)
c
	integer	i,j,k,ires
	integer	idx,ipocket
	integer	lout
	integer	nedge,ntetra,npockets
	integer	npoints,nvertex,ntrig
	integer	trig1,trig2,trig3,trig4
	integer	pair1,pair2,pair3,pair4,pair5,pair6
c
	integer redinfo(npointmax)
c
	integer edge(2,nedge_max)
	integer*1 edge_status(nedge_max)
	integer*2 edge_stat(nedge_max)
	integer*2 edge_stat2(nedge_max)
c
	integer trig_link(3,ntrig_max)
	integer	trig(3,ntrig_max)
c
	integer*1 trig_status(ntrig_max),trig_coef(ntrig_max)
	integer   tetra(4,ntetra_max),tetra_neighbour(4,ntetra_max)
	integer   tetra_link(4,ntetra_max)
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
	integer	tetra_pocket(ntetra_max)
c
	character name*4
	character resname*4,chname*1
	character label(natot)*54,tag(natot)*14
	character name_lc*4,res_lc*4,ch_lc*4
	character to_lower4*4,to_lower1
c
	common  /vertex_zone/   npoints,nvertex,redinfo
	common /xyz_vertex/     coord,radius,weight
	common /trig_zone/	ntrig,trig
	common /trig_stat/	trig_status,trig_coef
	common  /tetra_zone/    ntetra,tetra,tetra_neighbour
	common  /tetra_stat/    tetra_status,tetra_orient,
     1				tetra_nindex
	common  /edge_zone/	nedge,edge
	common  /edge_stat/	edge_status
	common  /links/         tetra_link,trig_link
	common  /pockets/       npockets,tetra_pocket
c
1	format(13x,a4,a4,a1,i4)
2	format('{ ',a14,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3,' { ',a14,' }',1x,a1,1x,f8.3,',',f8.3,',',
     2	f8.3)
3	format('{ ',a14,' }',1x,a1,1x,f8.3,',',f8.3,',',
     1	f8.3)
4	format('@vectorlist {Delcx} color=blue')
5	format('@vectorlist {AComplex1} color=red')
6	format('@vectorlist {Pck',i1,'} color=green')
7	format('@vectorlist {Pck',i2,'} color=green')
8	format('@vectorlist {Pck',i3,'} color=green')
9	format('@group {Geometry} ')
10	format(a4,a4,a1,1x,i4)
11	format('@subgroup {    Delaunay} dominant off')
12	format('@subgroup {    A-complex} dominant off')
13	format('@subgroup {    Pocket',i1,'} dominant off')
14	format('@subgroup {    Pocket',i2,'} dominant off')
15	format('@subgroup {    Pocket',i3,'} dominant off')
16	format('@vectorlist {AComplex2} color=cyan')
17	format('@vectorlist {AComplex3} color=magenta')
c
	write(lout,*) ' '
	write(lout,9)
	write(lout,*) ' '
	do 100 i = 1,nedge
		edge_stat(i) = 0
		edge_stat2(i) = 0
100	continue
c
c
	do 200 idx = 1,ntetra
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
		ipocket = tetra_pocket(idx)
		if(ipocket.eq.0) goto 150
c
		edge_stat(pair1) = ipocket
		edge_stat(pair2) = ipocket
		edge_stat(pair3) = ipocket
		edge_stat(pair4) = ipocket
		edge_stat(pair5) = ipocket
		edge_stat(pair6) = ipocket
c
150		continue
c
		if(tetra_status(idx).eq.1) then
			edge_stat2(pair1) = -1
			edge_stat2(pair2) = -1
			edge_stat2(pair3) = -1
			edge_stat2(pair4) = -1
			edge_stat2(pair5) = -1
			edge_stat2(pair6) = -1
		endif
c
200	continue
c
	do 250 idx = 1,ntrig
c
		if(trig_status(idx).eq.1) then
c
			pair1 = trig_link(1,idx)
			pair2 = trig_link(2,idx)
			pair3 = trig_link(3,idx)
c
			if(edge_stat2(pair1).eq.0) then
				edge_stat2(pair1) = -2
			endif
			if(edge_stat2(pair2).eq.0) then
				edge_stat2(pair2) = -2
			endif
			if(edge_stat2(pair3).eq.0) then
				edge_stat2(pair3) = -2
			endif
c
		endif
c
250	continue
c
	do 280 idx = 1,nedge
		if(edge_status(idx).eq.1.and.
     1			edge_stat2(idx).eq.0) then
			edge_stat2(idx) = -3
		endif
280	continue
c
	do 300 i = 1,npoints
		read(label(i),1) name,resname,chname,
     1		ires
		name_lc = to_lower4(name)
		ch_lc = to_lower1(chname)
		res_lc = to_lower4(resname)
		write(tag(i),10) name_lc,res_lc,ch_lc,
     1		ires
300	continue
c
c	Now write the Delaunay
c
	write(lout,*) ' '
	write(lout,11)
	write(lout,4)
c
	do 400 idx = 1,nedge
c
		i = edge(1,idx)
		j = edge(2,idx)
c
		write(lout,2) tag(i),'P',(coord(3*(i-1)+k),
     1		k=1,3),tag(i),'L',(coord(3*(i-1)+k),
     2		k=1,3)
		write(lout,3) tag(j),'L',(coord(3*(j-1)+k),
     1		k=1,3)
c
400	continue
c
c	write Alpha complex
c
c	We go through the edges three times : first if they
c	belong to tetrahedron, second if they belong to
c	singular triangles, and third if they are singular
c	themselves
c
	write(lout,*) ' '
	write(lout,12)
	write(lout,5)
	do 500 idx = 1,nedge
c
		if(edge_stat2(idx).ne.-1) goto 500
c
		i = edge(1,idx)
		j = edge(2,idx)
c
		write(lout,2) tag(i),'P',(coord(3*(i-1)+k),
     1		k=1,3),tag(i),'L',(coord(3*(i-1)+k),
     2		k=1,3)
		write(lout,3) tag(j),'L',(coord(3*(j-1)+k),
     1		k=1,3)
c
500	continue
	write(lout,16)
	do 600 idx = 1,nedge
c
		if(edge_stat2(idx).ne.-2) goto 600
c
		i = edge(1,idx)
		j = edge(2,idx)
c
		write(lout,2) tag(i),'P',(coord(3*(i-1)+k),
     1		k=1,3),tag(i),'L',(coord(3*(i-1)+k),
     2		k=1,3)
		write(lout,3) tag(j),'L',(coord(3*(j-1)+k),
     1		k=1,3)
c
600	continue
	write(lout,17)
	do 700 idx = 1,nedge
c
		if(edge_stat2(idx).ne.-3) goto 700
c
		i = edge(1,idx)
		j = edge(2,idx)
c
		write(lout,2) tag(i),'P',(coord(3*(i-1)+k),
     1		k=1,3),tag(i),'L',(coord(3*(i-1)+k),
     2		k=1,3)
		write(lout,3) tag(j),'L',(coord(3*(j-1)+k),
     1		k=1,3)
c
700	continue
c
c	Now write all the pockets
c
	do 900 ipocket = 1,min(10,npockets)
c
		write(lout,*) ' '
		if(ipocket.lt.10) then
			write(lout,13) ipocket
			write(lout,6) ipocket
		elseif(ipocket.lt.100) then
			write(lout,14) ipocket
			write(lout,7) ipocket
		else
			write(lout,15) ipocket
			write(lout,8) ipocket
		endif
		do 800 idx = 1,nedge
c
			if(edge_stat(idx).ne.ipocket) goto 800
c
			i = edge(1,idx)
			j = edge(2,idx)
c
			write(lout,2) tag(i),'P',(coord(3*(i-1)+k),
     1			k=1,3),tag(i),'L',(coord(3*(i-1)+k),
     2			k=1,3)
			write(lout,3) tag(j),'L',(coord(3*(j-1)+k),
     1			k=1,3)
c
800		continue
c
900	continue
c
	return
	end
c
c	kin_cpk.f		Version 1 5/14/2002	Patrice Koehl
c
c	This program generates a KineMage script to view a molecule in
c	CPK format
c
	subroutine kin_cpk(lout,label,radii,radius_h2o)
c
	include 'pocket.h'
c
	real*8	radius_h2o,radius
	real*8	coord(ncortot)
	real*8	radii(nresdef,20)
c
	integer	nat,i,j,k,natom,nseq,ires
	integer	lout
	integer	itype(nrestot),listatom(nrestot)
	integer	chain(nrestot),occup(natot)
c
	character atom(natot)*4,seq(nrestot)*4
	character name*4
	character resname*4,chname*1
	character label(natot)*54
	character name_lc*4,res_lc*4,ch_lc*4
	character to_lower4*4,to_lower1
c
	common /protein/ nseq,itype,chain,seq
	common /xyz/     coord,occup,natom,listatom
	common /names/   atom
c
1	format('@group {Balls} ')
2	format('@balllist {Ball_VdW} color=gold radius = ',f6.3)
3	format('@balllist {Ball_H2O} color=green radius = ',f6.3)
4	format('{ ',a4,a4,a1,1x,i4,' }',1x,f8.3,',',f8.3,',',
     1	f8.3)
5	format(21x,a1)
6	format('@subgroup {    CPK_VdW} dominant off')
7	format('@subgroup {    CPK_H2O} dominant off')
c
c	Start CPK model
c
	write(lout,1)
c
c	First define VdW model
c
	write(lout,6)
	nat = listatom(1)
c
	do 200 i = 2,nseq-1
		resname = seq(i)
		if(i.eq.1) resname = seq(i+1)
		if(i.eq.nseq) resname = seq(i-1)
		ires = i
		if(i.gt.1) ires = ires-1
		if(i.eq.nseq) ires = ires -1
		res_lc = to_lower4(resname)
		do 100 j = 1,listatom(i)
			nat = nat + 1
			if(atom(nat)(1:1).eq.'H') goto 100
			if(occup(nat).eq.0) goto 100
			if(j.eq.1) then
				read(label(nat),5) chname
				ch_lc = to_lower1(chname)
			endif
			name = atom(nat)
			name_lc = to_lower4(name)
			radius=radii(itype(i),j)
			write(lout,2) radius
			write(lout,4) name_lc,res_lc,ch_lc,ires,
     1				(coord(3*nat-3+k),k=1,3)
100		continue
200	continue
c
c	Now define "accessible" model
c
	write(lout,7)
c
	nat = listatom(1)
c
	do 400 i = 2,nseq-1
		resname = seq(i)
		if(i.eq.1) resname = seq(i+1)
		if(i.eq.nseq) resname = seq(i-1)
		ires = i
		if(i.gt.1) ires = ires-1
		if(i.eq.nseq) ires = ires -1
		res_lc = to_lower4(resname)
		do 300 j = 1,listatom(i)
			nat = nat + 1
			if(atom(nat)(1:1).eq.'H') goto 300
			if(occup(nat).eq.0) goto 300
			if(j.eq.1) then
				read(label(nat),5) chname
				ch_lc = to_lower1(chname)
			endif
			name = atom(nat)
			name_lc = to_lower4(name)
			radius=radii(itype(i),j)+radius_h2o
			write(lout,3) radius
			write(lout,4) name_lc,res_lc,ch_lc,ires,
     1				(coord(3*nat-3+k),k=1,3)
300		continue
400	continue
c
	return
	end
