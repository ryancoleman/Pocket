c	Readrad.f
c
c	This subroutine reads a parameter file containing all atomic
c	radii of all atoms of the 20 amino acids
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
	subroutine readrad(fname,radii,aspval,itest)
c
	include 'pocket.h'
c
	integer		i,iadd,ires,itest
	integer		residue
	integer		idihed(nresdef),iatom(nresdef)
	integer		atomtype(nresdef,20)
c
	real*8		radius,asp
	real*8		radii(nresdef,20)
	real*8		aspval(nresdef,20)
	real*8		atomcharge(nresdef,20)
c
	character	fname*30
	character	record*80,name*4,res*4
	character	nameatom(nresdef,20)*4
c
	common /default/ nameatom,atomtype,atomcharge,idihed,iatom
c
1	format(a)
2	format(a4,1x,a4,21x,f7.3,f7.3)
c
	itest = 1
c
	open(unit=1,file=fname,status='old')
c
100	read(1,1,end=200) record
	if(record(1:1).eq.'#'.or.record(1:1).eq.'r') goto 100
c
	if(record(41:41).ne.'.') then
		itest = 0
		read(record,2) res,name,radius
	else
		read(record,2) res,name,radius,asp
	endif
c
	ires = residue(res)
	if(ires.eq.0) goto 100
c
        iadd = 0
c
	if(name.eq.'H   ') name = 'HN'
	if(name.eq.'HT3') name = 'HN'
c
	if(name.eq.'HT1') then
		iadd = -1
		goto 40
	endif
	if(name.eq.'HT2') then
		iadd = 0
		goto 40
	endif
c
	if(name.eq.'OT1 ') name = 'O'
c
	if(name.eq.'OT2') then
		iadd = 0
		goto 40
	endif
c
	if(name.eq.'OXT') then
		iadd = 0
		goto 40
	endif
c
	if(ires.eq.4.and.name.eq.'CD') name = 'CD1'
	if(ires.eq.11.and.name.eq.'HG') name = 'HG1'
c
	do 20 i = 1,iatom(ires)
		if(name.eq.nameatom(ires,i)) then
			iadd = i
			goto 40
		endif
20	continue
c
40	continue
	if(iadd.eq.0.or.iadd.eq.-1) goto 100
	radii(ires,iadd) = radius
	if(itest.eq.1) aspval(ires,iadd) = asp
	goto 100
c
200	continue
c
	return
	end
