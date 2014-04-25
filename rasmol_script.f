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
C	Rasmol_script.f		Version 1 5/13/2002	Patrice Koehl
c
c	This subroutine generates a rasmol script for representing
c	the molecule studied, and its voids and pockets
c
	subroutine rasmol_script(fname,label,natom_pocket,
     1		atom_pocket)
c
	include 'pocket.h'
c
	integer	i,ipocket,istart,icolor
	integer	natom_pocket
c
	integer	atom_pocket(natot)
c
	character fname*50
	character text*30,atom*15,apocket*10
	character text_pocket*30,build_name*15,pocket_name*10
	character label(natot)*54
	character color_list(10)*9
c
1	format(a30)
2	format('select ',a15)
3	format('select selected or ',a15)
4	format('select ',a10)
5	format('color ',a9)
c
	data color_list/'blue','green','red','yellow','magenta',
     1	'cyan','orange','redorange','violet','purple'/
c
	open(unit=1,file=fname,status='unknown')
c
	ipocket = 0
	istart = 0
	do 100 i = 1,natom_pocket
c
		if(atom_pocket(i).lt.0) then
			if(ipocket.ne.0) then
				text=text_pocket(ipocket)
				write(1,1) text
			endif
			ipocket = -atom_pocket(i)
			istart = 1
			goto 100
		endif
c
		atom = build_name(label(atom_pocket(i)))
		if(istart.eq.1) then
			write(1,2) atom
			istart = 0
		else
			write(1,3) atom
		endif
c
100	continue
	text=text_pocket(ipocket)
	write(1,1) text
c
	write(1,*) 'restrict none'
	write(1,*) 'background white'
	write(1,*) 'select all'
	write(1,*) 'wireframe'
	write(1,*) 'backbone 50'
c
	do 200 i = 1,ipocket
c
		apocket = pocket_name(i)
		write(1,4) apocket
		icolor = mod(i,10)+1
		write(1,5) color_list(icolor)
		write(1,*) 'spacefill'
c
200	continue
c
	close(unit=1)
c
	return
	end
c
c	text_pocket.f		Version 1 5/13/2002	Patrice Koehl
c
c	This function generate the "define pocket" line for the
c	rasmol script
c
	function text_pocket(ipocket)
c
	integer	ipocket
c
	character	text_pocket*30
c
1	format('define pocket',i1,' selected')
2	format('define pocket',i2,' selected')
3	format('define pocket',i3,' selected')
4	format('define pocket',i4,' selected')
c
	if(ipocket.lt.10) then
		write(text_pocket,1) ipocket
	elseif(ipocket.lt.100) then
		write(text_pocket,2) ipocket
	elseif(ipocket.lt.1000) then
		write(text_pocket,3) ipocket
	elseif(ipocket.lt.10000) then
		write(text_pocket,4) ipocket
	endif
c
	return
	end
c
c	pocket_name.f		Version 1 5/13/2002	Patrice Koehl
c
c	This function generate the pocket name
c
	function pocket_name(ipocket)
c
	integer	ipocket
c
	character	pocket_name*10
c
1	format('pocket',i1)
2	format('pocket',i2)
3	format('pocket',i3)
4	format('pocket',i4)
c
	if(ipocket.lt.10) then
		write(pocket_name,1) ipocket
	elseif(ipocket.lt.100) then
		write(pocket_name,2) ipocket
	elseif(ipocket.lt.1000) then
		write(pocket_name,3) ipocket
	elseif(ipocket.lt.10000) then
		write(pocket_name,4) ipocket
	endif
c
	return
	end
c
c	Build_name.f	Version 1 5/13/2002	Patrice Koehl
c
c	This function generates the atom name from the PDb record
c	(where the atom name follows the rasmol format)
c
	function build_name(label)
c
	integer		i,j
c
	character*15 	build_name,name
	character	label*54
	character	atom*4,residue*4,chain*1,resnum*4
c
1	format(12x,a4,1x,a4,a1,a4)
c
	read(label,1) atom,residue,chain,resnum
c
	name=residue//resnum//chain//'.'//atom
c
	build_name='               '
c
	j = 0
	do 100 i = 1,14
		if(name(i:i).ne.' ') then
			j = j + 1
			build_name(j:j) = name(i:i)
		endif
100	continue
c
	return
	end
