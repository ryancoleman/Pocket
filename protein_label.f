c	Protein.f
c
c	This file contains a library of routines for dealing with
c	PDB files containing protein coordinates
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
c	Topology.f	Version 1 12/7/1995		Patrice Koehl
c
c	This file contains a blockdata which initialise the name of the 
c	standard amino acids that can be found in a protein, the number of
c	atoms and dihedrals which define each of these amino acids, as well
c	as the name, category and charge of each of these atoms.
c	Only the 20 standard amino acids are considered at this stage
c	Atom category are defined according to CHARMM19 : only polar
c	hydrogens are included
c
c	Atom types, according to CHARMM (used for VdW interactions) :
c
c	Hydrogens :
c			1.  H  (H which can H-bond to neutral atom)
c			2.  HC (H which can H-bond to charged atom)
c			3.  HA (aliphatic hydrogen)
c			4.  HT (TIPS3P water hydrogen)
c			5.  LP (ST2 lone pair)
c	Carbons
c			6.  CT   (aliphatic carbon)
c			7.  C    (carbonyl carbon)
c			8.  CH1E (extended carbon, with 1 hydrogen)
c			9.  CH2E (extended carbon, with 2 hydrogens)
c			10. CH3E (extended carbon, with 3 hydrogens)
c			11. CR1E (ring carbons)
c			12. CM   (C in carbonmonoxide)
c	Nitrogens
c			13. N	(peptide N with no hydrogen atoms attached)
c			14. NR  (N in aromatic ring with no hydrogen atoms)
c			15. NP	(Pyrole N)
c			16. NH1 (Peptide N bound to one hydrogen)
c			17. NH2 (Peptide N bound to two hydrogen atoms)
c			18. NH3 (Peptide N bound to three hydrogen atoms)
c			19. NC  (Guanidinium N bound to 2 hydrogens)
c
c	Oxygens
c			20. O   (carbonyl O)
c			21. OC  (carboxy O)
c			22. OH1 (hydroxy O)
c			23. OH2 (ST2 water oxygen)
c			24. OM  (O in carbonmonoxide) 
c			25. OT  (TIPS3P water oxygen)
c			26. OS  (esther oxygen)
c                       27. OXT (terminal oxygen at end of chain)
c
c	Sulfurs
c			28. S   (sulfur)
c			29. SH1E (extended atom S with one hydrogen)
c	
c	Iron
c			30. FE  (iron)
c
	subroutine topology_new
c
	include 'pocket.h'
c
	integer*1	exclusion(nresdef,400)
	integer*1	ex_i_iplus1(15,70),disu(7,7),disu2(7,7)
c
	integer		i,j,ntype
	integer		idihed(nresdef),iatom(nresdef)
	integer		atomtype(nresdef,20)
c
	real*8		atomcharge(nresdef,20)
c
	character	nameatom(nresdef,20)*4,nameres(nresdef)*4
	character	atom_name(natom_type)*4
c
	common /default/ nameatom,atomtype,atomcharge,idihed,iatom
	common /name/	 ntype,nameres
	common /exclude/ exclusion,ex_i_iplus1,disu,disu2
	common /atomtyp/ atom_name
c
	data (atom_name(i),i=1,30) /'H   ','HC  ','HA  ','HT  ','LP  ',
     1	'CT  ','C   ','CH1E','CH2E','CH3E','CR1E','CM  ','N   ','NR  ',
     2	'NP  ','NH1 ','NH2 ','NH3 ','NC2 ','O   ','OC  ','OH1 ','OH2 ',
     3  'OM  ','OT  ','OS  ','OXT ','S   ','SH1E','FE  '/
c
	data ntype/24/
c
	data (nameres(i),i=1,24)/'GLY','ALA','VAL','ILE','LEU','PHE',
     1		'PRO','MET','TRP','CYS','SER','THR','ASN','GLN','TYR',
     2		'HIS','ASP','GLU','LYS','ARG','HSD','HSP','NTER','CTER'/
c
	data (idihed(i),i=1,24)/0,0,1,2,2,2,1,3,2,1,2,2,3,4,2,2,2,3,
     1				5,7,2,2,0,0/
c
	data (iatom(i),i=1,24)/6,6,8,9,9,12,8,9,16,7,8,9,11,12,14,12,9,
     1				10,13,17,12,13,2,1/
c
c	1. Glycine : no sidechain, only 5 atoms to define the backbone
c	For sake of simplicity, a fictitious atoms is added in place of
c	CB, such that all residues (except NTER and CTER), have 6 atoms
c	in their backbone
c
	data (nameatom(1,i),i=1,6) /'N','HN','CA','X','C','O'/
	data (atomtype(1,i),i=1,6) /16,1,9,0,7,20/
	data (atomcharge(1,i),i=1,6)/-0.35,0.25,0.10,0,0.55,-0.55/
c
c	2. Alanine : only CB to define the sidechain
c
	data (nameatom(2,i),i=1,6) /'N','HN','CA','CB','C','O'/
	data (atomtype(2,i),i=1,6) /16,1,8,10,7,20/
	data (atomcharge(2,i),i=1,6)/-0.35,0.25,0.10,0.,0.55,-0.55/
c
c	3. Valine : 
c
	data (nameatom(3,i),i=1,8) /'N','HN','CA','CB','C','O',
     1		'CG1','CG2'/
	data (atomtype(3,i),i=1,8) /16,1,8,8,7,20,10,10/
	data (atomcharge(3,i),i=1,8) /-0.35,0.25,0.10,0,0.55,-0.55,
     1	0.,0./
c
c	4. Isoleucine :
c
	data (nameatom(4,i),i=1,9) /'N','HN','CA','CB','C','O','CG1',
     1	'CG2','CD1'/
	data (atomtype(4,i),i=1,9) /16,1,8,8,7,20,9,10,10/
	data (atomcharge(4,i),i=1,9) /-0.35,0.25,0.10,0.,0.55,-0.55,
     1	0.,0.,0./
c
c	5. Leucine :
c
	data (nameatom(5,i),i=1,9) /'N','HN','CA','CB','C','O','CG',
     1	'CD1','CD2'/
	data (atomtype(5,i),i=1,9) /16,1,8,9,7,20,8,10,10/
	data (atomcharge(5,i),i=1,9) /-0.35,0.25,0.10,0.,0.55,-0.55,
     1	0.,0.,0./
c
c	6. Phenylalanine :
c
	data (nameatom(6,i),i=1,12) /'N','HN','CA','CB','C','O','CG',
     1	'CD1','CD2','CE1','CE2','CZ'/
	data (atomtype(6,i),i=1,12) /16,1,8,9,7,20,7,11,11,11,11,11/
	data (atomcharge(6,i),i=1,12) /-0.35,0.25,0.10,0.,0.55,
     1	-0.55,0.,0.,0.,0.,0.,0./
c
c	7. Proline :
c
	data (nameatom(7,i),i=1,8) /'N','X','CA','CB','C','O','CG','CD'/
	data (atomtype(7,i),i=1,8) /13,0,8,9,7,20,9,9/
	data (atomcharge(7,i),i=1,8) /-0.2,0.,0.10,0.,0.55,-0.55,
     1	0.,0.10/
c
c	8. Methionine
c
	data (nameatom(8,i),i=1,9) /'N','HN','CA','CB','C','O','CG',
     1	'SD','CE'/
	data (atomtype(8,i),i=1,9) /16,1,8,9,7,20,9,27,10/
	data (atomcharge(8,i),i=1,9) /-0.35,0.25,0.10,0.,0.55,-0.55,
     1	0.06,-0.12,0.06/
c
c	9. Tryptophane
c
	data (nameatom(9,i),i=1,16) /'N','HN','CA','CB','C','O','CG',
     1	'CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2','HE1'/
	data (atomtype(9,i),i=1,16) /16,1,8,9,7,20,7,11,7,16,7,
     1	11,11,11,11,1/
	data (atomcharge(9,i),i=1,16) /-0.35,0.25,0.10,0.,0.55,-0.55,
     1	-0.03,0.06,0.10,-0.36,-0.04,-0.03,0.,0.,0.,0.3/
c
c	10. Cysteine
c
	data (nameatom(10,i),i=1,7) /'N','HN','CA','CB','C','O','SG'/
	data (atomtype(10,i),i=1,7) /16,1,8,9,7,20,28/
	data (atomcharge(10,i),i=1,7) /-0.35,0.25,0.10,0.19,0.55,
     1	-0.55,-0.19/
c
c	11. Serine
c
	data (nameatom(11,i),i=1,8) /'N','HN','CA','CB','C','O','OG',
     1	'HG1'/
	data (atomtype(11,i),i=1,8) /16,1,8,9,7,20,22,1/
	data (atomcharge(11,i),i=1,8) /-0.35,0.25,0.10,0.25,0.55,-0.55,
     1	-0.65,0.40/
c
c	12. Threonin
c
	data (nameatom(12,i),i=1,9) /'N','HN','CA','CB','C','O','CG2',
     1	'OG1','HG1'/
	data (atomtype(12,i),i=1,9) /16,1,8,8,7,20,10,22,1/
	data (atomcharge(12,i),i=1,9) /-0.35,0.25,0.10,0.25,0.55,-0.55,
     1	0.,-0.65,0.40/
c
c	13. Asparagine
c
	data (nameatom(13,i),i=1,11) /'N','HN','CA','CB','C','O','CG',
     1	'OD1','ND2','HD21','HD22'/
	data (atomtype(13,i),i=1,11) /16,1,8,9,7,20,7,20,17,1,1/
	data (atomcharge(13,i),i=1,11) /-0.35,0.25,0.10,0,0.55,-0.55,
     1	0.55,-0.55,-0.60,0.30,0.30/
c
c	14. Glutamine
c
	data (nameatom(14,i),i=1,12) /'N','HN','CA','CB','C','O','CG',
     1	'CD','OE1','NE2','HE21','HE22'/
	data (atomtype(14,i),i=1,12) /16,1,8,9,7,20,9,7,20,17,1,1/
	data (atomcharge(14,i),i=1,12) /-0.35,0.25,0.10,0.,0.55,-0.55,
     1	0.,0.55,-0.55,-0.6,0.3,0.3/
c
c	15. Tyrosine
c
	data (nameatom(15,i),i=1,14) /'N','HN','CA','CB','C','O','CG',
     1	'CD1','CD2','CE1','CE2','CZ','OH','HH'/
	data (atomtype(15,i),i=1,14) /16,1,8,9,7,20,7,11,11,11,11,7,
     1	22,1/
	data (atomcharge(15,i),i=1,14)/-0.35,0.25,0.10,0.,0.55,-0.55,
     1	0.,0.,0.,0.,0.,0.25,-0.65,0.40/
c
c	16. Histidine (ND1 protonated)
c
	data (nameatom(16,i),i=1,12) /'N','HN','CA','CB','C','O','CG',
     1	'ND1','CD2','CE1','NE2','HD1'/
	data (atomtype(16,i),i=1,12) /16,1,8,9,7,20,7,16,11,11,14,1/
	data (atomcharge(16,i),i=1,12) /-0.35,0.25,0.10,0.,0.55,-0.55,
     1	0.10,-0.40,0.10,0.30,-0.40,0.30/
c
c	17. Aspartic Acid
c
	data (nameatom(17,i),i=1,9) /'N','HN','CA','CB','C','O','CG',
     1	'OD1','OD2'/
	data (atomtype(17,i),i=1,9) /16,1,8,9,7,20,7,21,21/
	data (atomcharge(17,i),i=1,9) /-0.35,0.25,0.10,-0.16,0.55,-0.55,
     1	0.36,-0.6,-0.6/
c
c	18. Glutamic Acid
c
	data (nameatom(18,i),i=1,10) /'N','HN','CA','CB','C','O','CG',
     1	'CD','OE1','OE2'/
	data (atomtype(18,i),i=1,10) /16,1,8,9,7,20,9,7,21,21/
	data (atomcharge(18,i),i=1,10) /-0.35,0.25,0.10,0.,0.55,-0.55,
     1	-0.16,0.36,-0.6,-0.6/
c
c	19. Lysine
c
	data (nameatom(19,i),i=1,13) /'N','HN','CA','CB','C','O','CG',
     1  'CD','CE','NZ','HZ1','HZ2','HZ3'/
	data (atomtype(19,i),i=1,13) /16,1,8,9,7,20,9,9,9,18,2,2,2/
	data (atomcharge(19,i),i=1,13) /-0.35,0.25,0.10,0.,0.55,-0.55,
     1	0.,0.,0.25,-0.30,0.35,0.35,0.35/
c
c	20. Arginine
c
	data (nameatom(20,i),i=1,17) /'N','HN','CA','CB','C','O','CG',
     1	'CD','NE','CZ','HE','NH1','NH2','HH11','HH12','HH21','HH22'/
	data (atomtype(20,i),i=1,17) /16,1,8,9,7,20,9,9,16,7,1,16,16,
     1	2,2,2,2/
	data (atomcharge(20,i),i=1,17) /-0.35,0.25,0.10,0.,0.55,
     1	-0.55,0.,0.10,-0.4,0.50,0.30,-0.45,-0.45,0.35,0.35,0.35,0.35/
c
c	21. Histidine (HSD : NE2 protonated)
c
	data (nameatom(21,i),i=1,12) /'N','HN','CA','CB','C','O','CG',
     1	'ND1','CD2','CE1','NE2','HE2'/
	data (atomtype(21,i),i=1,12) /16,1,8,9,7,20,7,14,11,11,16,1/
	data (atomcharge(21,i),i=1,12) /-0.35,0.25,0.10,0.,0.55,-0.55,
     1	0.10,-0.40,0.10,0.30,-0.40,0.30/
c
c	22. Histidine (HSC : protonated ND1 and NE2)
c
	data (nameatom(22,i),i=1,13) /'N','HN','CA','CB','C','O','CG',
     1	'ND1','CD2','CE1','NE2','HD1','HE2'/
	data (atomtype(22,i),i=1,13) /16,1,8,9,7,20,7,14,11,11,16,1,1/
	data (atomcharge(22,i),i=1,13) /-0.35,0.25,0.10,0.10,0.55,-0.55,
     1	0.15,-0.30,0.20,0.45,-0.30,0.35,0.35/
c
c	23. Nter
c
	data (nameatom(23,i),i=1,2) /'HT1','HT2'/
	data (atomtype(23,i),i=1,2) /2,2/
	data (atomcharge(23,i),i=1,2) /0.35,0.35/
c
c	24. Cter
c
	data (nameatom(24,i),i=1,1) /'OXT'/
	data (atomtype(24,i),i=1,1) /21/
	data (atomcharge(24,i),i=1,1) /-0.57/
c
c	Definition of contingency table for each type of amino acid :
c	Possible non bonded interactions between two atoms i and j
c	are classified according to :
c	0 : 1-2 (bond) 
c	1 : 1-3 (angle) connectivities 
c		(in both cases, the two atoms are 
c	    	excluded from the list of non bonded interactions)
c	2 : 1-4 interactions (i.e 1 torsion angle) : can be reduced both
c	    for VdW (different LJ potentials) as well as electrostatics
c	    (E14FAC)
c	3 : all other interactions
c	-1 : refers to i, or j (or both) being non existant
c
c	Values are provided only for i<j
c
c	1. Gly :
c			      X(4)    O(6)          1 2 3 4 5 6
c			      |	      |		  1 x 0 0-1 1 2
c			N(1)--CA(3)--C(5)         2 x x 1-1 2 3
c			|                         3 x x x-1 0 1
c			HN(2)                     4 x x x x-1-1
c						  5 x x x x x 0
c
	data (exclusion(1,i),i=1,15) /0,0,-1,1,2,1,-1,2,3,-1,0,1,
     1	-1,-1,0/
c
c	2. Ala :					1 2 3 4 5 6
c			         CB(4)  O(6)	      1 x 0 0 1 1 2
c				 |      ||            2 x x 1 2 2 3
c			   N(1)--CA(3)--C(5)          3 x x x 0 0 1
c			   |                          4 x x x x 1 2
c			   HN(2)		      5 x x x x x 0
c
	data (exclusion(2,i),i=1,15) /0,0,1,1,2,1,2,2,3,0,0,1,1,2,0/
c
c	3. Val				1 2 3 4 5 6 7 8
c		CG1(7)  CG2(8)        1 x 0 0 1 1 2 2 2
c		   \    /             2 x x 1 2 2 3 3 3
c		    CB(4)  O(6)       3 x x x 0 0 1 1 1
c		     |     ||         4 x x x x 1 2 0 0
c	      N(1)--CA(3)--C(5)       5 x x x x x 0 2 2
c	      |			      6 x x x x x x 3 3
c	      HN(2) 		      7 x x x x x x x 1
c
	data (exclusion(3,i),i=1,28) /0,0,1,1,2,2,2,1,2,2,3,3,3,0,0,
     1	1,1,1,1,2,0,0,0,2,2,3,3,1/
c
c	4. Ile				1 2 3 4 5 6 7 8 9
c		CD1(9)                1 x 0 0 1 1 2 2 2 3
c		|                     2 x x 1 2 2 3 3 3 3
c		CG1(7)  CG2(8)        3 x x x 0 0 1 1 1 2
c		  \    /              4 x x x x 1 2 0 0 1
c		    CB(4)  O(6)       5 x x x x x 0 2 2 3
c		     |     ||         6 x x x x x x 3 3 3
c	      N(1)--CA(3)--C(5)       7 x x x x x x x 1 0
c	      |			      8 x x x x x x x x 2
c	      HN(2) 
c
	data (exclusion(4,i),i=1,36) /0,0,1,1,2,2,2,3,1,2,2,3,3,3,3,
     1	0,0,1,1,1,2,1,2,0,0,1,0,2,2,3,3,3,3,1,0,2/
c
c	5. Leu				1 2 3 4 5 6 7 8 9
c		CD1(8)  CD2(9)        1 x 0 0 1 1 2 2 3 3
c		   \    /             2 x x 1 2 2 3 3 3 3
c		    CG(7)             3 x x x 0 0 1 1 2 2
c		     |                4 x x x x 1 2 0 1 1
c		    CB(4)  O(6)       5 x x x x x 0 2 3 3
c		     |     ||         6 x x x x x x 3 3 3
c	      N(1)--CA(3)--C(5)       7 x x x x x x x 0 0
c	      |			      8 x x x x x x x x 1
c	      HN(2) 
c
	data (exclusion(5,i),i=1,36) /0,0,1,1,2,2,3,3,1,2,2,3,3,3,3,
     1	0,0,1,1,2,2,1,2,0,1,1,0,2,3,3,3,3,3,0,0,1/
c
c	6. Phe                                            1 1 1
c					1 2 3 4 5 6 7 8 9 0 1 2   
c		     CZ (12)          1 x 0 0 1 1 2 2 3 3 3 3 3
c		   /    \             2 x x 1 2 2 3 3 3 3 3 3 3
c	      	CE1(10) CE2(11)       3 x x x 0 0 1 1 2 2 3 3 3
c		|         |           4 x x x x 1 2 0 1 1 2 2 3
c		CD1(8)  CD2(9)        5 x x x x x 0 2 3 3 3 3 3
c		   \    /             6 x x x x x x 3 3 3 3 3 3
c		    CG(7)             7 x x x x x x x 0 0 1 1 1
c		     |                8 x x x x x x x x 1 0 1 1
c		    CB(4)  O(6)       9 x x x x x x x x x 1 0 1
c		     |     ||        10 x x x x x x x x x x 1 0
c	      N(1)--CA(3)--C(5)      11 x x x x x x x x x x x 0
c	      |			     
c	      HN(2) 
c
	data (exclusion(6,i),i=1,66) /0,0,1,1,2,2,3,3,3,3,3,1,2,2,3,3,3,
     1  3,3,3,3,0,0,1,1,2,2,3,3,3,1,2,0,1,1,2,2,3,0,2,3,3,3,3,3,3,3,3,3,
     2  3,3,0,0,1,1,1,1,0,1,1,1,0,1,1,0,0/
c
c	7. Pro                           1 2 3 4 5 6 7 8
c		 CG(7)                 1 x-1 0 1 1 2 1 0
c	       /    \                  2 x x-1-1-1-1-1-1
c	      CD(8)  CB(4)  O(6)       3 x x x 0 0 1 1 1 
c             |	     |     ||          4 x x x x 1 2 0 1
c	      N(1)--CA(3)--C(5)        5 x x x x x 0 2 2
c	      |			       6 x x x x x x 3 3
c             X(2)                     7 x x x x x x x 0
c
	data (exclusion(7,i),i=1,28) /-1,0,1,1,2,1,0,6*-1,0,0,1,1,
     1	1,1,2,0,1,0,2,2,3,3,0/
c
c	8. Met
c		    CE(9)
c		     |                  1 2 3 4 5 6 7 8 9
c		    SD(8)             1 x 0 0 1 1 2 2 3 3
c                    |                2 x x 1 2 2 3 3 3 3
c		    CG(7)             3 x x x 0 0 1 1 2 3
c		     |                4 x x x x 1 2 0 1 2
c		    CB(4)  O(6)       5 x x x x x 0 2 3 3
c		     |     ||         6 x x x x x x 3 3 3
c	      N(1)--CA(3)--C(5)       7 x x x x x x x 0 1
c	      |			      8 x x x x x x x x 0
c	      HN(2) 
c
	data (exclusion(8,i),i=1,36) / 0,0,1,1,2,2,3,3,1,2,2,3,3,3,3,
     1  0,0,1,1,2,3,1,2,0,1,2,0,2,3,3,3,3,3,0,1,0/
c
c	9. Trp                                                  1 1 1 1 1 1 1
c                                             1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6
c                           CZ2(13)         1 x 0 0 1 1 2 2 3 3 3 3 3 3 3 3 3
c                         /     \           2 x x 1 2 2 3 3 3 3 3 3 3 3 3 3 3
c     HE1(16)--NE1(10)--CE2(11)  CH2(15)    3 x x x 0 0 1 1 2 2 3 3 3 3 3 3 3
c               |        |       |          4 x x x x 1 2 0 1 1 2 2 2 3 3 3 3
c              CD1(8)   CD2(9)   CZ3(14)    5 x x x x x 0 2 3 3 3 3 3 3 3 3 3
c                \     /  \     /           6 x x x x x x 3 3 3 3 3 3 3 3 3 3
c                  CG(7)    CE3(12)         7 x x x x x x x 0 0 1 1 1 1 1 1 2
c	            |                       8 x x x x x x x x 1 0 1 1 1 1 1 1
c		   CB(4)  O(6)              9 x x x x x x x x x 1 0 0 1 1 1 2
c		    |     ||               10 x x x x x x x x x x 0 1 1 1 1 0
c	     N(1)--CA(3)--C(5)             11 x x x x x x x x x x x 1 0 1 1 1
c	     |			           12 x x x x x x x x x x x x 1 0 1 3
c	     HN(2)                         13 x x x x x x x x x x x x x 1 0 2
c					   14 x x x x x x x x x x x x x x 0 3
c					   15 x x x x x x x x x x x x x x x 3
c
	data (exclusion(9,i),i=1,120) /0,0,1,1,2,2,3,3,3,3,3,3,3,
     1	3,3,1,2,2,3,3,3,3,3,3,3,3,3,3,3,0,0,1,1,2,2,3,3,3,3,3,3,
     2	3,1,2,0,1,1,2,2,2,3,3,3,3,0,2,3,3,3,3,3,3,3,3,3,10*3,0,
     3  0,6*1,2,1,0,6*1,1,0,0,1,1,1,2,0,1,1,1,1,0,1,0,1,1,1,1,0,
     4	1,3,1,0,2,0,3,3/
c
c	10. Cys
c                                       1 2 3 4 5 6 7
c		    SG(7)             1 x 0 0 1 1 2 2 
c		     |                2 x x 1 2 2 3 3 
c		    CB(4)  O(6)       3 x x x 0 0 1 1 
c		     |     ||         4 x x x x 1 2 0
c	      N(1)--CA(3)--C(5)       5 x x x x x 0 2 
c	      |			      6 x x x x x x 3
c	      HN(2) 		     
c
	data (exclusion(10,i),i=1,21) /0,0,1,1,2,2,1,2,2,3,3,0,0,1,1,
     1	1,2,0,0,2,3/
c
c	11. Ser
c		    HG1(8)              1 2 3 4 5 6 7 8
c		     |                1 x 0 0 1 1 2 2 3
c		    OG(7)             2 x x 1 2 2 3 3 3
c		     |                3 x x x 0 0 1 1 2
c		    CB(4)  O(6)       4 x x x x 1 2 0 1
c		     |     ||         5 x x x x x 0 2 3
c	      N(1)--CA(3)--C(5)       6 x x x x x x 3 3
c	      |			      7 x x x x x x x 0
c	      HN(2) 		     
c
	data (exclusion(11,i),i=1,28) /0,0,1,1,2,2,3,1,2,2,3,3,3,0,0,
     1	1,1,2,1,2,0,1,0,2,3,3,3,0/
c
c	12. Thr
c					1 2 3 4 5 6 7 8 9
c		HG1(9)                1 x 0 0 1 1 2 2 2 3
c		|                     2 x x 1 2 2 3 3 3 3
c		OG1(8)  CG2(7)        3 x x x 0 0 1 1 1 2
c		  \    /              4 x x x x 1 2 0 0 1
c		    CB(4)  O(6)       5 x x x x x 0 2 2 3
c		     |     ||         6 x x x x x x 3 3 3
c	      N(1)--CA(3)--C(5)       7 x x x x x x x 1 2
c	      |			      8 x x x x x x x x 0
c	      HN(2) 
c
	data (exclusion(12,i),i=1,36) /0,0,1,1,2,2,2,3,1,2,2,3,3,3,3,
     1	0,0,1,1,1,2,1,2,0,0,1,0,2,2,3,3,3,3,1,2,0/
c
c	13. Asparagine
c                                                         1 1
c		   HD21(10) HD22(11)    1 2 3 4 5 6 7 8 9 0 1
c	               \  /	      1 x 0 0 1 1 2 2 3 3 3 3
c		OD1(8)  ND2(9)        2 x x 1 2 2 3 3 3 3 3 3
c		   \    /             3 x x x 0 0 1 1 2 2 3 3
c		    CG(7)             4 x x x x 1 2 0 1 1 2 2
c		     |                5 x x x x x 0 2 3 3 3 3
c		    CB(4)  O(6)       6 x x x x x x 3 3 3 3 3
c		     |     ||         7 x x x x x x x 0 0 1 1
c	      N(1)--CA(3)--C(5)       8 x x x x x x x x 1 2 2
c	      |			      9 x x x x x x x x x 0 0
c	      HN(2)                  10 x x x x x x x x x x 1
c
	data (exclusion(13,i),i=1,55) /0,0,1,1,2,2,3,3,3,3,1,2,2,3,3,
     1	3,3,3,3,0,0,1,1,2,2,3,3,1,2,0,1,1,2,2,0,2,3,3,3,3,3,3,3,3,3,0,
     2	0,1,1,1,2,2,0,0,1/
c
c	14. Gln
c
c                                                         1 1 1 
c		   HE21(11) HE22(12)    1 2 3 4 5 6 7 8 9 0 1 2
c	               \  /	      1 x 0 0 1 1 2 2 3 3 3 3 3
c		OE1(9)  NE2(10)       2 x x 1 2 2 3 3 3 3 3 3 3
c		   \    /             3 x x x 0 0 1 1 2 3 3 3 3
c		    CD(8)             4 x x x x 1 2 0 1 2 2 3 3
c                    |                5 x x x x x 0 2 3 3 3 3 3
c		    CG(7)             6 x x x x x x 3 3 3 3 3 3
c		     |                7 x x x x x x x 0 1 1 2 2
c		    CB(4)  O(6)       8 x x x x x x x x 0 0 1 1
c		     |     ||         9 x x x x x x x x x 1 2 2
c	      N(1)--CA(3)--C(5)      10 x x x x x x x x x x 0 0
c	      |			     11 x x x x x x x x x x x 1
c	      HN(2)                  
c
	data (exclusion(14,i),i=1,66) /0,0,1,1,2,2,3,3,3,3,3,1,2,2,
     1	3,3,3,3,3,3,3,0,0,1,1,2,3,3,3,3,1,2,0,1,2,2,3,3,0,2,3,3,3,3,
     2  3,3,3,3,3,3,3,0,1,1,2,2,0,0,1,1,1,2,2,0,0,1/
c
c	15. Tyr
c
c		     HH (14)
c		     |
c	             OH (13)                              1 1 1 1 1
c	             |			1 2 3 4 5 6 7 8 9 0 1 2 3 4 
c		     CZ (12)          1 x 0 0 1 1 2 2 3 3 3 3 3 3 3
c		   /    \             2 x x 1 2 2 3 3 3 3 3 3 3 3 3 
c	      	CE1(10) CE2(11)       3 x x x 0 0 1 1 2 2 3 3 3 3 3 
c		|         |           4 x x x x 1 2 0 1 1 2 2 3 3 3 
c		CD1(8)  CD2(9)        5 x x x x x 0 2 3 3 3 3 3 3 3
c		   \    /             6 x x x x x x 3 3 3 3 3 3 3 3
c		    CG(7)             7 x x x x x x x 0 0 1 1 1 3 3
c		     |                8 x x x x x x x x 1 0 1 1 2 3
c		    CB(4)  O(6)       9 x x x x x x x x x 1 0 1 2 3
c		     |     ||        10 x x x x x x x x x x 1 0 1 2
c	      N(1)--CA(3)--C(5)      11 x x x x x x x x x x x 0 1 2
c	      |			     12 x x x x x x x x x x x x 0 1
c	      HN(2)                  13 x x x x x x x x x x x x x 0
c
	data (exclusion(15,i),i=1,91) /0,0,1,1,2,2,7*3,1,2,2,9*3,0,
     1	0,1,1,2,2,5*3,1,2,0,1,1,2,2,3,3,3,0,2,7*3,8*3,0,0,3*1,3,3,
     2  1,0,1,1,2,3,1,0,1,2,3,1,0,1,2,0,1,2,0,1,0/
c
c	16. Histidine
c							    1 1 1
c	      NE2(11)-- CE1(10)           1 2 3 4 5 6 7 8 9 0 1 2 
c	       |         |              1 x 0 0 1 1 2 2 3 3 3 3 3 
c	      CD2(9)    ND1(8)--HD1(12) 2 x x 1 2 2 3 3 3 3 3 3 3 
c	         \     /                3 x x x 0 0 1 1 2 2 3 3 3  
c		  CG(7)                 4 x x x x 1 2 0 1 1 2 2 2
c		     |                  5 x x x x x 0 2 3 3 3 3 3
c		    CB(4)  O(6)         6 x x x x x x 3 3 3 3 3 3  
c		     |     ||           7 x x x x x x x 0 0 1 1 1  
c	      N(1)--CA(3)--C(5)         8 x x x x x x x x 1 0 1 0  
c	      |			        9 x x x x x x x x x 1 0 2  
c	      HN(2)                    10 x x x x x x x x x x 0 1  
c                                      11 x x x x x x x x x x x 2
c
	data (exclusion(16,i),i=1,66) /0,0,1,1,2,2,5*3,1,2,2,7*3,0,0,
     1  1,1,2,2,3,3,3,1,2,0,1,1,2,2,2,0,2,11*3,0,0,1,1,1,1,0,1,0,1,0,
     2	2,0,1,2/
c
c	17. Asp
c                                       1 2 3 4 5 6 7 8 9
c		OD1(8)  OD2(9)        1 x 0 0 1 1 2 2 3 3    
c		   \    /             2 x x 1 2 2 3 3 3 3 
c		    CG(7)             3 x x x 0 0 1 1 2 2    
c		     |                4 x x x x 1 2 0 1 1 
c		    CB(4)  O(6)       5 x x x x x 0 2 3 3 
c		     |     ||         6 x x x x x x 3 3 3     
c	      N(1)--CA(3)--C(5)       7 x x x x x x x 0 0    
c	      |			      8 x x x x x x x x 1    
c	      HN(2)                  
c
	data (exclusion(17,i),i=1,36) /0,0,1,1,2,2,3,3,1,2,2,4*3,
     1	0,0,1,1,2,2,1,2,0,1,1,0,2,3,3,3,3,3,0,0,1/
c
c	18. Glu                                           1
c                                       1 2 3 4 5 6 7 8 9 0
c				      1 x 0 0 1 1 2 2 3 3 3
c		OE1(9)  OE2(10)       2 x x 1 2 2 3 3 3 3 3 
c		   \    /             3 x x x 0 0 1 1 2 3 3 
c		    CD(8)             4 x x x x 1 2 0 1 2 2 
c                    |                5 x x x x x 0 2 3 3 3 
c		    CG(7)             6 x x x x x x 3 3 3 3 
c		     |                7 x x x x x x x 0 1 1 
c		    CB(4)  O(6)       8 x x x x x x x x 0 0 
c		     |     ||         9 x x x x x x x x x 1 
c	      N(1)--CA(3)--C(5)      
c	      |			     
c	      HN(2)                  
c
	data (exclusion(18,i),i=1,45) /0,0,1,1,2,2,3,3,3,1,2,2,5*3,
     1	0,0,1,1,2,3,3,1,2,0,1,2,2,0,2,3,3,3,3,3,3,3,0,1,1,0,0,1/
c
c	19. Lys
c                                                          1 1 1 1
c	   HZ1(11)  HZ2(12) HZ3(13)      1 2 3 4 5 6 7 8 9 0 1 2 3
c	          \   |    /           1 x 0 0 1 1 2 2 3 3 3 3 3 3
c		    NZ(10)             2 x x 1 2 2 3 3 3 3 3 3 3 3
c		     |                 3 x x x 0 0 1 1 2 3 3 3 3 3
c                   CE(9)              4 x x x x 1 2 0 1 2 3 3 3 3
c		     |                 5 x x x x x 0 2 3 3 3 3 3 3
c		    CD(8)              6 x x x x x x 3 3 3 3 3 3 3
c                    |                 7 x x x x x x x 0 1 2 3 3 3
c		    CG(7)              8 x x x x x x x x 0 1 2 2 2
c		     |                 9 x x x x x x x x x 0 1 1 1
c		    CB(4)  O(6)       10 x x x x x x x x x x 0 0 0
c		     |     ||         11 x x x x x x x x x x x 1 1
c	      N(1)--CA(3)--C(5)       12 x x x x x x x x x x x x 1
c	      |			     
c	      HN(2)                  
c
	data (exclusion(19,i),i=1,78) /0,0,1,1,2,2,6*3,1,2,2,8*3,0,0,1,
     1  1,2,5*3,1,2,0,1,2,4*3,0,2,6*3,7*3,0,1,2,3,3,3,0,1,2,2,2,0,1,1,
     2  1,0,0,0,1,1,1/
c
c	20. Arg
c                                                          1 1 1 1 1 1 1 1
c   HH11(14) HH12(15) HH21(16) HH22(17)  1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7
c	\    /            \    /       1 x 0 0 1 1 2 2 3 3 3 3 3 3 3 3 3 3
c         NH1(12)         NH2(13)      2 x x 1 2 2 3 3 3 3 3 3 3 3 3 3 3 3
c	       \        /              3 x x x 0 0 1 1 2 3 3 3 3 3 3 3 3 3
c                 CZ(10)               4 x x x x 1 2 0 1 2 3 3 3 3 3 3 3 3
c		     |                 5 x x x x x 0 2 3 3 3 3 3 3 3 3 3 3
c		    NE(9)--HE(11)      6 x x x x x x 3 3 3 3 3 3 3 3 3 3 3
c                    |                 7 x x x x x x x 0 1 2 2 3 3 3 3 3 3
c		    CD(8)              8 x x x x x x x x 0 1 1 2 2 3 3 3 3
c                    |                 9 x x x x x x x x x 0 0 1 1 2 2 2 2
c		    CG(7)             10 x x x x x x x x x x 1 0 0 1 1 1 1
c		     |                11 x x x x x x x x x x x 2 2 3 3 3 3
c		    CB(4)  O(6)       12 x x x x x x x x x x x x 1 0 0 2 2
c		     |     ||         13 x x x x x x x x x x x x x 2 2 0 0
c	      N(1)--CA(3)--C(5)       14 x x x x x x x x x x x x x x 1 3 3
c	      |			      15 x x x x x x x x x x x x x x x 3 3
c	      HN(2)                   16 x x x x x x x x x x x x x x x x 1
c
	data (exclusion(20,i),i=1,136) /0,0,1,1,2,2,10*3,1,2,2,12*3,
     1	0,0,1,1,2,9*3,1,2,0,1,2,8*3,0,2,10*3,11*3,0,1,2,2,6*3,0,1,1,2,2,
     1  3,3,3,3,0,0,1,1,2,2,2,2,1,0,0,1,1,1,1,2,2,3,3,3,3,
     1	1,0,0,2,2,2,2,0,0,1,3,3,3,3,1/
c
c	21. HSD (NE2 protonated)
c							    1 1 1
c    HE2(12)--NE2(11)-- CE1(10)           1 2 3 4 5 6 7 8 9 0 1 2 
c	       |         |              1 x 0 0 1 1 2 2 3 3 3 3 3 
c	      CD2(9)    ND1(8)          2 x x 1 2 2 3 3 3 3 3 3 3 
c	         \     /                3 x x x 0 0 1 1 2 2 3 3 3  
c		  CG(7)                 4 x x x x 1 2 0 1 1 2 2 3
c		     |                  5 x x x x x 0 2 3 3 3 3 3
c		    CB(4)  O(6)         6 x x x x x x 3 3 3 3 3 3  
c		     |     ||           7 x x x x x x x 0 0 1 1 2  
c	      N(1)--CA(3)--C(5)         8 x x x x x x x x 1 0 1 2  
c	      |			        9 x x x x x x x x x 1 0 1  
c	      HN(2)                    10 x x x x x x x x x x 0 1  
c                                      11 x x x x x x x x x x x 0
c
	data (exclusion(21,i),i=1,66) /0,0,1,1,2,2,3,3,3,3,3,1,2,2,
     1  7*3,0,0,1,1,2,2,3,3,3,1,2,0,1,1,2,2,3,0,2,3,3,3,3,3,3,3,3,
     1  3,3,3,0,0,1,1,2,1,0,1,2,1,0,1,0,1,0/
c
c	22. HSP (both HD1 and NE2 protonated)
c							    1 1 1 1
c    HE2(13)--NE2(11)-- CE1(10)           1 2 3 4 5 6 7 8 9 0 1 2 3
c	       |         |              1 x 0 0 1 1 2 2 3 3 3 3 3 3
c	      CD2(9)    ND1(8)--HD1(12) 2 x x 1 2 2 3 3 3 3 3 3 3 3
c	         \     /                3 x x x 0 0 1 1 2 2 3 3 3 3
c		  CG(7)                 4 x x x x 1 2 0 1 1 2 2 2 3
c		     |                  5 x x x x x 0 2 3 3 3 3 3 3
c		    CB(4)  O(6)         6 x x x x x x 3 3 3 3 3 3 3
c		     |     ||           7 x x x x x x x 0 0 1 1 1 2
c	      N(1)--CA(3)--C(5)         8 x x x x x x x x 1 0 1 0 2
c	      |			        9 x x x x x x x x x 1 0 2 1
c	      HN(2)                    10 x x x x x x x x x x 0 1 1
c                                      11 x x x x x x x x x x x 2 0
c                                      12 x x x x x x x x x x x x 3
c
	data (exclusion(22,i),i=1,78) /0,0,1,1,2,2,6*3,1,2,2,8*3,0,0,
     1	1,1,2,2,4*3,1,2,0,1,1,2,2,2,3,0,2,13*3,0,0,1,1,1,2,1,0,1,0,2,
     2  1,0,2,1,0,1,1,2,0,3/
c
c	23. Nter :
c
c	HT1                            1 2
c	   \                         1 x 1
c	   /
c	HT2
c
	data (exclusion(23,i),i=1,1)/1/
c
c	24. Cter : no table required, since this pseudo residue only contains
c	one atom
c
c	Now exclusion list for the backbones of two consecutive residues :
c
c	Case 1 : Glycine - Glycine :
c
c             X(4)   O(6)         X(4')   O(6')       1' 2' 3' 4' 5' 6'
c             |      ||           |       ||        1 2  3  3 -1  3  3
c	N(1)--CA(3)--C(5)--N(1')--CA(3')--C(5')     2 3  3  3 -1  3  3
c	|                  |                        3 1  2  2 -1  3  3 
c	HN(2)              HN(2')                   4-1 -1 -1 -1 -1 -1
c                                                   5 0  1  1 -1  2  3
c						    6 1  2  2 -1  3  3
c
	data (ex_i_iplus1(1,i),i=1,36)/2,3,3,-1,3,3,3,3,3,-1,3,3,1,2,2,
     1	-1,3,3,-1,-1,-1,-1,-1,-1,0,1,1,-1,2,3,1,2,2,-1,3,3/
c
c	Case 2 : Glycine - Proline :
c
c				CG(7')                1' 2' 3' 4' 5' 6' 7' 8'
c                            /     \                1 2 -1  3  3  3  3  3  3
c             X(4)   O(6) CD(8') CB(4')   O(6')     2 3 -1  3  3  3  3  3  3
c	      |	     ||    |        /      ||       3 1 -1  2  3  3  3  3  2
c	N(1)--CA(3)--C(5)--N(1')--CA(3')--C(5')     4-1 -1 -1 -1 -1 -1 -1 -1
c	|                  |                        5 0 -1  1  2  2  3  2  1
c	HN(2)              X(2')                    6 1 -1  2  3  3  3  3  2
c
	data (ex_i_iplus1(2,i),i=1,48)/2,-1,6*3,3,-1,6*3,1,-1,2,4*3,2,
     1	8*-1,0,-1,1,2,2,3,2,1,1,-1,2,4*3,2/
c
c	Case 3 : Glycine - OT2 :
c
c	      X(4)   OT1(6)                               1'
c	      |	     ||                                1  2
c	N(1)--CA(3)--C(5)--OT2(1')                     2  3
c	|                                              3  1
c	HN(2)                                          4 -1 
c                                                      5  0
c						       6  1
c
	data (ex_i_iplus1(3,i),i=1,6) /2,3,1,-1,0,1/
c
c	Case 4 : Glycine - Non Gly, non Pro amino acids
c
c	      X(4)   O(6)         CB(4')  O(6')       1' 2' 3' 4' 5' 6'
c	      |	     ||           |       ||        1 2  3  3  3  3  3
c	N(1)--CA(3)--C(5)--N(1')--CA(3')--C(5')     2 3  3  3  3  3  3
c	|                  |                        3 1  2  2  3  3  3 
c	HN(2)              HN(2')                   4-1 -1 -1 -1 -1 -1
c						    5 0  1  1  2  2  3
c                                                   6 1  2  2  3  3  3
c
	data (ex_i_iplus1(4,i),i=1,36) /2,5*3,6*3,1,2,2,3*3,6*-1,0,
     1	1,1,2,2,3,1,2,2,3,3,3/
c
c	Case 5 : NTER - Gly residue
c
c	HT1(1)      X(4')   OT1(6')                  1' 2' 3' 4' 5' 6'
c          \        |       ||                     1 0  1  1 -1  2  3
c	     N(1')--CA(3')--C(5')                  2 0  1  1 -1  2  3
c	   / |                                          
c    HT2(2)  HT3(2')                        
c
	data (ex_i_iplus1(5,i),i=1,12) /0,1,1,-1,2,3,0,1,1,-1,2,3/
c
c	Case 6 : Nter - Other residues
c
c	HT1(1)      CB(4')  OT1(6')                  1' 2' 3' 4' 5' 6'
c          \        |       ||                     1 0  1  1  2  2  3
c	     N(1')--CA(3')--C(5')                  2 0  1  1  2  2  3
c	   / |                                          
c    HT2(2)  HT3(2')                        
c
	data (ex_i_iplus1(6,i),i=1,12) /0,1,1,2,2,3,0,1,1,2,2,3/
c
c	Case 7 : Non Gly - Gly 
c
c             CB(4)  O(6)         X(4')   O(6')       1' 2' 3' 4' 5' 6'
c	      |	     ||           |       ||        1 2  3  3 -1  3  3
c	N(1)--CA(3)--C(5)--N(1')--CA(3')--C(5')     2 3  3  3 -1  3  3
c	|                  |                        3 1  2  2 -1  3  3
c	HN(2)              HN(2')                   4 2  3  3 -1  3  3
c                                                   5 0  1  1 -1  2  3
c                                                   6 1  2  2 -1  3  3
c
	data (ex_i_iplus1(7,i),i=1,36) /2,3,3,-1,3,3,3,3,3,-1,3,3,1,2,2,
     1	-1,3,3,2,3,3,-1,3,3,0,1,1,-1,2,3,1,2,2,-1,3,3/
c
c	Case 8 : Non Gly - Pro
c
c				CG(7')                1' 2' 3' 4' 5' 6' 7' 8'
c                            /     \                1 2 -1  3  3  3  3  3  3
c	      CB(4)  O(6) CD(8') CB(4')   O(6')     2 3 -1  3  3  3  3  3  3
c             |      ||    |        /      ||       3 1 -1  2  3  3  3  3  2
c	N(1)--CA(3)--C(5)--N(1')--CA(3')--C(5')     4 2 -1  3  3  3  3  3  3
c	|                  |                        5 0 -1  1  2  2  3  2  1
c	HN(2)              X(2')                    6 1 -1  2  3  3  3  3  2
c
	data (ex_i_iplus1(8,i),i=1,48) /2,-1,6*3,3,-1,6*3,1,-1,2,3,3,3,
     1	3,2,2,-1,6*3,0,-1,1,2,2,3,2,1,1,-1,2,3,3,3,3,2/
c
c	Case 9 : Non Gly - OT2
c                                                         1'
c	      CB(4)  OT1(6)                            1  2
c             |      ||                                2  3
c	N(1)--CA(3)--C(5)--OT2(1')                     3  1
c	|                                              4  2
c	HN(2)                                          5  0
c                                                      6  1
c
	data (ex_i_iplus1(9,i),i=1,6) /2,3,1,2,0,1/
c
c	Case 10 : Non Gly - Non Gly, non Pro amino acids
c
c             CB(4)  O(6)         CB(4')  O(6')       1' 2' 3' 4' 5' 6'
c             |      ||           |       ||        1 2  3  3  3  3  3
c	N(1)--CA(3)--C(5)--N(1')--CA(3')--C(5')     2 3  3  3  3  3  3
c	|                  |                        3 1  2  2  3  3  3 
c	HN(2)              HN(2')                   4 2  3  3  3  3  3
c                                                   5 0  1  1  2  2  3
c						    6 1  2  2  3  3  3
c
	data (ex_i_iplus1(10,i),i=1,36) /2,5*3,6*3,1,2,2,3*3,2,5*3,
     1	0,1,1,2,2,3,1,2,2,3,3,3/
c
c	Case 11 : Pro - Gly 
c
c         CG(7)
c        /   \
c      CD(8) CB(4)  O(6)         X(4')   O(6')       1' 2' 3' 4' 5' 6'
c	|     |	     ||           |       ||        1 2  3  3 -1  3  3
c	N(1)--CA(3)--C(5)--N(1')--CA(3')--C(5')     2-1 -1 -1 -1 -1 -1
c	|                  |                        3 1  2  2 -1  3  3
c	X(2)              HN(2')                    4 2  3  3 -1  3  3
c                                                   5 0  1  1 -1  2  3
c                                                   6 1  2  2 -1  3  3
c                                                   7 3  3  3 -1  3  3
c                                                   8 3  3  3 -1  3  3
c
	data (ex_i_iplus1(11,i),i=1,48) /2,2*3,-1,2*3,6*-1,1,2,2,-1,
     1	2*3,2,2*3,-1,2*3,0,1,1,-1,2,3,1,2,2,-1,3,3,3,3,3,-1,3,3,3,3,3,
     2  -1,3,3/
c
c	Case 12 : Pro - Pro 
c
c         CG(7)               CG(7')
c        /   \              /    \
c      CD(8) CB(4)  O(6)   CD(8') CB(4')   O(6')      1' 2' 3' 4' 5' 6' 7' 8'
c	|     |	     ||    |      |       ||        1 2 -1  3  3  3  3  3  3
c	N(1)--CA(3)--C(5)--N(1')--CA(3')--C(5')     2-1 -1 -1 -1 -1 -1 -1 -1
c	|                  |                        3 1 -1  2  3  3  3  3  2
c	X(2)               X(2')                    4 2 -1  3  3  3  3  3  3
c                                                   5 0 -1  1  2  2  3  2  1
c                                                   6 1 -1  2  3  3  3  3  2
c                                                   7 3 -1  3  3  3  3  3  3   
c                                                   8 3 -1  3  3  3  3  3  3
c
	data (ex_i_iplus1(12,i),i=1,64) /2,-1,6*3,8*-1,1,-1,2,3,3,3,
     1	3,2,2,-1,6*3,0,-1,1,2,2,3,2,1,1,-1,2,3,3,3,3,2,3,-1,3,3,3,3,3,3,
     2  3,-1,3,3,4*3/
c
c	Case 13 : Pro - Non Pro, Non Gly
c
c         CG(7)
c        /   \
c      CD(8) CB(4)  O(6)         CB(4')   O(6')       1' 2' 3' 4' 5' 6'
c	|     |	     ||           |       ||        1 2  3  3  3  3  3
c	N(1)--CA(3)--C(5)--N(1')--CA(3')--C(5')     2-1 -1 -1 -1 -1 -1
c	|                  |                        3 1  2  2  3  3  3
c	X(2)              HN(2')                    4 2  3  3  3  3  3
c                                                   5 0  1  1  2  2  3
c                                                   6 1  2  2  3  3  3
c                                                   7 3  3  3  3  3  3
c                                                   8 3  3  3  3  3  3
c
	data (ex_i_iplus1(13,i),i=1,48) /2,5*3,6*-1,1,2,2,3*3,2,
     1	5*3,0,1,1,2,2,3,1,2,2,3,3,3,12*3/
c
c	Case 14 : Pro - OT2 
c
c         CG(7)
c        /   \
c      CD(8) CB(4)  O(6)               1' 
c	|     |	     ||              1 2  
c	N(1)--CA(3)--C(5)--OXT(1')   2-1 
c	|                            3 1 
c	X(2)                         4 2 
c                                    5 0 
c                                    6 1 
c                                    7 3 
c                                    8 3
c
	data (ex_i_iplus1(14,i),i=1,8) /2,-1,1,2,0,1,3,3/
c
c	Finally, special case for disulphide bridges
c
c       HN(2')
c       |
c	N(1')--CA(3')--C(5')
c	       |       ||
c	       CB(4')  O(6')              1' 2' 3' 4' 5' 6' 7'
c              |                       1  3  3  3  3  3  3  3
c              SG(7')                  2  3  3  3  3  3  3  3
c              |                       3  3  3  3  3  3  3  2
c	       SG(7)                   4  3  3  3  2  3  3  1
c	       |                       5  3  3  3  3  3  3  3
c	       CB(4)  O(6)             6  3  3  3  3  3  3  3
c	       |     ||                7  3  3  2  1  3  3  0
c	N(1)--CA(3)--C(5)    
c	| 		    
c	HN(2) 		     
c
	data ((disu(i,j),j=1,7),i=1,7)/7*3,7*3,6*3,2,3*2,2,3,3,1,7*3,
     1	7*3,3,3,2,1,3,3,0/
c
c	             SG(7)--SG(7')                    1' 2' 3' 4' 5' 6' 7'
c	           /              \                 1 2  3  3  3  3  3  3
c	       CB(4)  O(6)        CB(4')  O(6')     2 3  3  3  3  3  3  3
c	       |     ||           |       ||        3 1  2  2  3  3  3  2
c	N(1)--CA(3)--C(5)--N(1')--CA(3')--C(5')     4 2  3  3  2  3  3  1
c	| 		   |                        5 0  1  1  2  3  3  3
c	HN(2) 		   HN(2')                   6 1  2  2  3  3  3  3
c                                                   7 3  3  2  1  3  3  0
c
	data ((disu2(i,j),j=1,7),i=1,7)/2,6*2,7*3,1,2,2,3,3,3,2,2,3,3,
     1	2,3,3,1,0,1,1,2,3,3,3,1,2,2,3,3,3,3,3,3,2,1,3,3,0/ 
c
	end
c	Defchain.f		Version 1 5/5/1992	Patrice Koehl
c
c	This program counts the number of chain in the protein
c
	subroutine defchain(fname,nchain,chname)
c
	integer	nchain
	integer	idum2,nres,ntest
c
	character fname*30,record*80,chname(10)*1,ch1
c
1	format(a)
2       format(22x,i4)
c
	open(unit=1,file=fname,status='old')
c
	ch1 	= 'X'
	nchain 	= 0
	nres	= 0
	ntest 	= -10
c
100     read(1,1,end=200) record
	if(record(1:6).eq.'ENDMDL') goto 200
        if(record(1:6).ne.'ATOM  ') goto 100
	if(record(27:27).ne.' ') goto 100
	if(record(22:22).ne.ch1) then
		nchain 			= nchain + 1
		chname(nchain) 		= record(22:22)
		ch1 			= chname(nchain)
	endif
	read(record,2) idum2
	if(idum2.ne.ntest) then
c                if(ntest.ne.-10.and.idum2-ntest.ne.1) then
c			write(6,*) 'Warning : non consecutive residues !!'
c			write(6,*) 'Residues : ',ntest,idum2
c                endif
                ntest = idum2
		nres  = nres +1
	endif
	goto 100
200	continue
c
	close(unit = 1)
c
	return
	end
c	extractseq.f		Version 1 24/7/1995	Patrice Koehl
c
c	this subroutine scans the pdb file to get the sequence of the
c	protein
c	Two 'pseudo' residues are added :
c	Nter and Cter (for HT1 and HT2, and for OT2, respectively)
c
	subroutine extractseq(fname)
c
	include 'pocket.h'
c
	integer	nseq,i,itest,residue
	integer	itype(nrestot),chain(nrestot)
c
	character fname*50
	character resname*4,record*80
	character seq(nrestot)*4
c
	common /protein/ nseq,itype,chain,seq
c
1	format(17x,a4,1x,i4)
2	format(a)
c
	open(unit=1,file=fname,status='old')
c
	nseq  = 1
	seq(1) = 'NTER'
	itest = -10
100     read(1,2,end=200) record
	if(record(1:6).eq.'ENDMDL') goto 200
        if(record(1:6).ne.'ATOM  ') goto 100
        read(record,1) resname,i
	if(i.ne.itest) then
		nseq = nseq + 1
		seq(nseq) = resname
		itest = i
	endif
	goto 100
200	continue
	nseq = nseq + 1
	seq(nseq) = 'CTER'
c
	close(unit=1)
c
        do 300 i = 1,nseq
                itype(i) = residue(seq(i))
300     continue
c
	return
	end
c	Prepare.f		Version 1 21/7/1995	Patrice Koehl
c
c	This subroutine prepares the molecule : 
c		- gives name of the atoms
c		- gives number of dihedral angles, number of atoms
c		  (total, and per residue)
c		- build up exclusion lists
c
	subroutine prepare
c
c	pocket.h contains the dimension of the arrays
c
	include 'pocket.h'
c
	integer		i,j,nseq,ndihed,natom,nco
	integer		itype(nrestot),listdihed(nrestot)
	integer		idihed(nresdef),iatom(nresdef)
	integer		atomtype(nresdef,20),listatom(nrestot)
	integer		ljtype(natot),ires(natot)
	integer		chain(nrestot),occup(natot)
c
	real*8		atomcharge(nresdef,20)
	real*8		dihed(nangtot),coord(ncortot)
	real*8		charge(natot)
c
	character 	atom(natot)*4,seq(nrestot)*4
	character	nameatom(nresdef,20)*4
c
	common /protein/ nseq,itype,chain,seq
	common /angle/   dihed,ndihed,listdihed
	common /xyz/	 coord,occup,natom,listatom
	common /param_ene/ charge,ljtype,ires
	common /names/	 atom
	common /default/ nameatom,atomtype,atomcharge,idihed,iatom
c
	ndihed = 0
	natom  = 0
c
	do 200 i = 1,nseq
c
c		1. Store all info for each residue :
c
		do 100 j = 1,iatom(itype(i))
			atom(natom+j)   = nameatom(itype(i),j)
			charge(natom+j) = atomcharge(itype(i),j)
			ljtype(natom+j) = atomtype(itype(i),j)
100		continue
c
c		2. Special case Nter and Cter residues :
c
c		a) For first residue, change charge of N and H, as 
c		well as name of H
c
		if(i.eq.2) then
			atom(natom+2)   = 'HT3'
			charge(natom+1) = -0.3d0
			charge(natom+2) = 0.35d0
			charge(natom+3) = 0.25d0
			ljtype(natom+2) = 2
		endif
c
c		b) for last residue, change charge of C, as well as
c		charge, name and type of O :
c
		if(i.eq.nseq-1) then
			nco = 5
			atom(natom+nco+1)   = 'OT1'
			charge(natom+nco)   = 0.14d0
			charge(natom+nco+1) = -0.57d0
			ljtype(natom+nco+1) = 12
		endif
c
		listatom(i)  = iatom(itype(i))
		if(i.gt.1.and.i.lt.nseq) then
			listdihed(i) = 3+idihed(itype(i))
		else
			listdihed(i) = idihed(itype(i))
		endif
c
		ndihed       = ndihed + listdihed(i)
		natom        = natom + listatom(i)
c
200	continue
c
c	write(6,*) ' '
c	write(6,*) ' Number of residue              : ',nseq-2
c	write(6,*) ' Number of atoms (with polar H) : ',natom
c	write(6,*) ' '
c
	return
	end
c	Residue.f		Version 1 24/7/1992	Patrice koehl
c
c       This functions gives a number (from 1 to ntype) to each amino acid
c
        integer function residue(aa)
c
	include 'pocket.h'
c
        character*4     aa
        character*4     nameres(nresdef)
c
        integer i,ntype
c
        common /name/   ntype,nameres
c
	if(aa.eq.'NTER') then
		residue = 23
		return
	endif
c
	if(aa.eq.'CTER') then
		residue = 24
		return
	endif
c
        do 100 i = 1,ntype
                if(aa.eq.nameres(i)) then
                        residue = i
                        return
                endif
100     continue
c
	residue = 0
c
        return
	end
c	Atompos.f		Version 1 24/7/1995	Patrice Koehl
c
c	This subroutine, knowing the name and the residue number of an
c	atom, gives its position in the array of atoms build for the 
c	calculation
c	Input :
c		- ires	: residu which contains the atom
c		- name	: name of the atom
c	Output :
c		- ipos	: atom number
c
	subroutine atompos(ires,name,ipos,iadd)
c
	include 'pocket.h'
c
c	pocket.h contains the dimension of the arrays
c
	integer		i,j,nseq,natom,ires,ipos,iadd
	integer		itype(nrestot),listatom(nrestot)
	integer		iatom(nresdef),idihed(nresdef)
	integer		atomtype(nresdef,20)
	integer		chain(nrestot),occup(natot)
c
	real*8		coord(ncortot)
	real*8		atomcharge(nresdef,20)
c
	character 	name*4,seq(nrestot)*4
	character	nameatom(nresdef,20)*4
c
	common /protein/ nseq,itype,chain,seq
	common /xyz/     coord,occup,natom,listatom
	common /default/ nameatom,atomtype,atomcharge,idihed,iatom
c
	ipos = 0
	do 10 i = 1,ires - 1
		ipos = ipos + listatom(i)
10	continue
c
        iadd = 0
c
	if(name.eq.'H   ') name = 'HN'
c
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
		iadd = 1 + listatom(ires)
		goto 40
	endif
c
	if(name.eq.'OXT') then
		iadd = 1 + listatom(ires)
		goto 40
	endif
c
	if(itype(ires).eq.4.and.name.eq.'CD') name = 'CD1'
c
	if(itype(ires).eq.11.and.name.eq.'HG') name = 'HG1'
c
	j = itype(ires)
	do 20 i = 1,iatom(j)
		if(name.eq.nameatom(j,i)) then
			iadd = i
			goto 40
		endif
20	continue
c
40	continue
c
	ipos = ipos + iadd
c
	return
	end
c	Addh.f		Version 1 28/8/1995		Patrice Koehl
c
c	This subroutine adds the missing polar hydrogens on sidechains
c
	subroutine addh(ityp,fulres,nfulres,checkat)
c
	integer	ityp,nfulres
	integer	checkat(nfulres)
c
	real*8	fulres(3*nfulres)
c
	if(ityp.eq.9) then
		call trph(fulres,nfulres,checkat)
	elseif(ityp.eq.11) then
		call serh(fulres,nfulres,checkat)
	elseif(ityp.eq.12) then
		call thrh(fulres,nfulres,checkat)
	elseif(ityp.eq.13) then
		call asnh(fulres,nfulres,checkat)
	elseif(ityp.eq.14) then
		call glnh(fulres,nfulres,checkat)
	elseif(ityp.eq.15) then
		call tyrh(fulres,nfulres,checkat)
	elseif(ityp.eq.16) then
		call hish(fulres,nfulres,checkat)
	elseif(ityp.eq.19) then
		call lysh(fulres,nfulres,checkat)
	elseif(ityp.eq.20) then
		call argh(fulres,nfulres,checkat)
	elseif(ityp.eq.21) then
		call hsdh(fulres,nfulres,checkat)
	elseif(ityp.eq.22) then
		call hsph(fulres,nfulres,checkat)
	endif
c
	return
	end
c	Buildh.f	Version 1 18/8/1995		Patrice Koehl
c
c	This set of subroutines adds missing polar hydrogens on some sidechains
c
c	1. Build HE1 on tryptophan
c
	subroutine trph(fulres,nfulres,checkat)
c
	integer	i,nfulres
	integer	checkat(nfulres)
c
	real*8	ang,dist,pif
	real*8	pa(3),pb(3),pc(3),pd(3)
	real*8	fulres(3*nfulres)
c
	pif = acos(-1.d0)/180.d0
c
	if(checkat(16).eq.1) return
	do 100 i = 1,3
		pa(i) = fulres(21+i)
		pb(i) = fulres(27+i)
		pc(i) = fulres(30+i)
100	continue
c
	dist = 1.0d0
	ang  = 126.d0*pif
	call sp2(pa,pb,pc,pd,ang,dist)
c
	do 200 i = 1,3
		fulres(45+i) = pd(i)
200	continue
c
	return
	end
c
c	2. Build HG1 on serine
c
	subroutine serh(fulres,nfulres,checkat)
c
	integer	i,nfulres
	integer	checkat(nfulres)
c
	real*8	ang,dist,tor,pif
	real*8	pa(3),pb(3),pc(3),pd(3)
	real*8	fulres(3*nfulres)
c
	pif = acos(-1.d0)/180.d0
c
	if(checkat(8).eq.1) return
c
	do 100 i = 1,3
		pa(i) = fulres(6+i)
		pb(i) = fulres(9+i)
		pc(i) = fulres(18+i)
100	continue
c
	dist = 1.0d0
	ang  = 110.d0*pif
	tor = 180.d0*pif
	call buildfour(pa,pb,pc,pd,tor,ang,dist)
c
	do 200 i = 1,3
		fulres(21+i) = pd(i)
200	continue
c
	return
	end
c
c	3. Build HG1 on threonin
c
	subroutine thrh(fulres,nfulres,checkat)
c
	integer	i,nfulres
	integer	checkat(nfulres)
c
	real*8	ang,dist,tor,pif
	real*8	pa(3),pb(3),pc(3),pd(3)
	real*8	fulres(3*nfulres)
c
	pif = acos(-1.d0)/180.d0
c
	if(checkat(9).eq.1) return
c
	do 100 i = 1,3
		pa(i) = fulres(6+i)
		pb(i) = fulres(9+i)
		pc(i) = fulres(21+i)
100	continue
c
	dist = 1.0d0
	ang  = 110.d0*pif
	tor  = 180.d0*pif
	call buildfour(pa,pb,pc,pd,tor,ang,dist)
c
	do 200 i = 1,3
		fulres(24+i) = pd(i)
200	continue
c
	return
	end
c
c	4. Build HD21 and HD22 on asparagin
c
	subroutine asnh(fulres,nfulres,checkat)
c
	integer	i,nfulres
	integer	checkat(nfulres)
c
	real*8	ang,dist,tor,pif
	real*8	pa(3),pb(3),pc(3),pd(3),pe(3)
	real*8	fulres(3*nfulres)
c
	pif = acos(-1.d0)/180.d0
c
	do 100 i = 1,3
		pa(i) = fulres(9+i)
		pb(i) = fulres(18+i)
		pc(i) = fulres(24+i)
100	continue
c
	dist = 1.0d0
	ang  = 120.d0*pif
	tor  = 180.d0*pif
	call buildfour(pa,pb,pc,pd,tor,ang,dist)
c
	call sp2(pb,pc,pd,pe,ang,dist)
c
	if(checkat(10).eq.0) then
		do 200 i = 1,3
			fulres(27+i) = pd(i)
200		continue
	endif
	if(checkat(11).eq.0) then
		do 300 i = 1,3
			fulres(30+i) = pe(i)
300		continue
	endif
c
	return
	end
c
c	5. Build HE21 and HE22 on glutamin
c
	subroutine glnh(fulres,nfulres,checkat)
c
	integer	i,nfulres
	integer	checkat(nfulres)
c
	real*8	ang,dist,tor,pif
	real*8	pa(3),pb(3),pc(3),pd(3),pe(3)
	real*8	fulres(3*nfulres)
c
	pif = acos(-1.d0)/180.d0
c
	do 100 i = 1,3
		pa(i) = fulres(18+i)
		pb(i) = fulres(21+i)
		pc(i) = fulres(27+i)
100	continue
c
	dist = 1.0d0
	ang  = 120.d0*pif
	tor  = 180.d0*pif
	call buildfour(pa,pb,pc,pd,tor,ang,dist)
c
	call sp2(pb,pc,pd,pe,ang,dist)
c
	if(checkat(11).eq.0) then
		do 200 i = 1,3
			fulres(30+i) = pd(i)
200		continue
	endif
	if(checkat(12).eq.0) then
		do 300 i = 1,3
			fulres(33+i) = pe(i)
300		continue
	endif
c
	return
	end
c
c	6. Build H of hydroxyl group of tyrosine
c
	subroutine tyrh(fulres,nfulres,checkat)
c
	integer	i,nfulres
	integer	checkat(nfulres)
c
	real*8	ang,dist,tor,pif
	real*8	pa(3),pb(3),pc(3),pd(3)
	real*8	fulres(3*nfulres)
c
	pif = acos(-1.d0)/180.d0
c
	if(checkat(14).eq.1) return
c
	do 100 i = 1,3
		pa(i) = fulres(27+i)
		pb(i) = fulres(33+i)
		pc(i) = fulres(36+i)
100	continue
c
	dist = 1.09d0
	tor  = 0.d0
	ang  = 120.d0*pif
	call buildfour(pa,pb,pc,pd,tor,ang,dist)
c
	do 200 i = 1,3
		fulres(39+i) = pd(i)
200	continue
c
	return
	end
c
c	7. Build HD1 on histidine
c
	subroutine hish(fulres,nfulres,checkat)
c
	integer	i,nfulres
	integer	checkat(nfulres)
c
	real*8	ang,dist,pif
	real*8	pa(3),pb(3),pc(3),pd(3)
	real*8	fulres(3*nfulres)
c
	pif = acos(-1.d0)/180.d0
c
	if(checkat(12).eq.1) return
c
	do 100 i = 1,3
		pa(i) = fulres(18+i)
		pb(i) = fulres(21+i)
		pc(i) = fulres(27+i)
100	continue
c
	dist = 1.0d0
	ang  = 125.d0*pif
	call sp2(pa,pb,pc,pd,ang,dist)
c
	do 200 i = 1,3
		fulres(33+i) = pd(i)
200	continue
c
	return
	end
c
c	8. Add three protons on lysine
c
	subroutine lysh(fulres,nfulres,checkat)
c
	integer	i,nfulres,ip
	integer	checkat(nfulres)
c
	real*8	ang,dist,pif,tor
	real*8	pa(3),pb(3),pc(3),pd(3),pe(3),pf(3)
	real*8	fulres(3*nfulres)
c
	pif = acos(-1.d0)/180.d0
c
	do 100 i = 1,3
		pa(i) = fulres(21+i)
		pb(i) = fulres(24+i)
		pc(i) = fulres(27+i)
100	continue
c
	dist = 1.01d0
	tor  = 180.d0*pif
	ang  = 112.d0*pif
	call buildfour(pa,pb,pc,pd,tor,ang,dist)
c
	ip = 1
	call nsp3(pb,pc,pd,pe,ang,ang,ip,dist)
	ip = 2
	call nsp3(pb,pc,pd,pf,ang,ang,ip,dist)
c
	if(checkat(11).eq.0.or.checkat(12).eq.0.or.
     1	checkat(13).eq.0) then
		do 200 i = 1,3
			fulres(30+i) = pd(i)
			fulres(33+i) = pe(i)
			fulres(36+i) = pf(i)
200		continue
	endif
c
	return
	end
c
c	9. Add protons on all terminal N on arginine sidechain
c
	subroutine argh(fulres,nfulres,checkat)
c
	integer i,nfulres
	integer	checkat(nfulres)
c
	real*8	ang,dist,pif,tor
	real*8	pa(3),pb(3),pc(3),pd(3),pe(3),pf(3),pg(3),ph(3)
	real*8	fulres(3*nfulres)
c
	pif = acos(-1.d0)/180.d0
c
	do 100 i = 1,3
		pa(i) = fulres(21+i)
		pb(i) = fulres(24+i)
		pc(i) = fulres(27+i)
100	continue
c
	dist = 1.01d0
	ang  = 123.d0*pif
	call sp2(pa,pb,pc,pd,ang,dist)
c
	do 200 i = 1,3
		pa(i) = pb(i)
		pb(i) = pc(i)
		pc(i) = fulres(33+i)
200	continue
c
	dist = 1.0d0
	ang  = 120.d0*pif
	tor  = 0
	call buildfour(pa,pb,pc,pe,tor,ang,dist)
	call sp2(pb,pc,pe,pf,ang,dist)
c
	do 300 i = 1,3
		pc(i) = fulres(36+i)
300	continue
c
	call buildfour(pa,pb,pc,pg,tor,ang,dist)
	call sp2(pb,pc,pg,ph,ang,dist)
c
	if(checkat(11).eq.0) then
		do 400 i = 1,3
			fulres(30+i) = pd(i)
400		continue
	endif
	if(checkat(14).eq.0.or.checkat(15).eq.0) then
		do 500 i = 1,3
			fulres(39+i) = pe(i)
			fulres(42+i) = pf(i)
500		continue
	endif
	if(checkat(16).eq.0.or.checkat(17).eq.0) then
		do 600 i = 1,3
			fulres(45+i) = pg(i)
			fulres(48+i) = ph(i)
600		continue
	endif
c
	return
	end
c
c	10. Add HE2 on HSD
c
	subroutine hsdh(fulres,nfulres,checkat)
c
	integer	i,nfulres
	integer	checkat(nfulres)
c
	real*8	ang,dist,pif
	real*8	pa(3),pb(3),pc(3),pd(3)
	real*8	fulres(3*nfulres)
c
	pif = acos(-1.d0)/180.d0
c
	if(checkat(12).ne.0) return
c
	do 100 i = 1,3
		pa(i) = fulres(24+i)
		pb(i) = fulres(30+i)
		pc(i) = fulres(27+i)
100	continue
c
	dist = 1.0d0
	ang  = 125.d0*pif
	call sp2(pa,pb,pc,pd,ang,dist)
c
	do 200 i = 1,3
		fulres(33+i) = pd(i)
200	continue
c
	return
	end
c
c	11. Add HD1 and HE2 on HSP
c
	subroutine hsph(fulres,nfulres,checkat)
c
	integer	i,nfulres
	integer	checkat(nfulres)
c
	real*8	ang,dist,pif
	real*8	pa(3),pb(3),pc(3),pd(3),pe(3)
	real*8	fulres(3*nfulres)
c
	pif = acos(-1.d0)/180.d0
c
	do 100 i = 1,3
		pa(i) = fulres(18+i)
		pb(i) = fulres(21+i)
		pc(i) = fulres(27+i)
100	continue
c
	dist = 1.0d0
	ang  = 125.d0*pif
	call sp2(pa,pb,pc,pd,ang,dist)
c
	do 200 i = 1,3
		pa(i) = fulres(24+i)
		pb(i) = fulres(30+i)
200	continue
c
	call sp2(pa,pb,pc,pe,ang,dist)
c
	if(checkat(12).eq.0) then
		do 300 i = 1,3
			fulres(33+i) = pd(i)
300		continue
	endif
	if(checkat(13).eq.0) then
		do 400 i = 1,3
			fulres(36+i) = pe(i)
400		continue
	endif
c
	return
	end
c	nsp3.f		version 1 8/3/1990		Patrice Koehl
c
c	This subroutine adds a proton H4 on a carbon with a sp3
c	conformation
c
c				X
c				!
c			C1 ---- C2 ---- C3
c				!
c				H4
c
c	Input of the program :
c				ang1	: angle C1-C2-H4 
c				ang2	: angle C2-C3-H4 
c				ip	: 1 for H4, 2 for X
c				dist	: length of C2-H4 (or C2-X)
c				p1	: C1
c				p2	: C2
c				p3	: C3
c	Output of the program :
c				p4	: H4 (or X)
c
	subroutine nsp3(p1,p2,p3,p4,ang1,ang2,ip,dist)
c
	real*8	p1(3),p2(3),p3(3),p4(3),ang1,ang2,dist
	real*8	v1(3),v2(3),v3(3),b(3),u1(3),dd(3),det
	real*8	dist1,dist2,d21(3),d23(3),pn(3),pm(3)
	real*8	alpha
c
	integer i,ip
c
	call diffvect(p2,p1,d21)
	call diffvect(p2,p3,d23)
	call normvect(d21,dist1)
	call normvect(d23,dist2)
	call crossvect(d21,d23,u1)
c
	do 10 i = 1,3
		v1(i) = dist*cos(ang1)*d21(i)/dist1
		v2(i) = dist*cos(ang2)*d23(i)/dist2
10	continue
c
	call addvect(v1,p2,pn)
	call addvect(v2,p2,pm)
c
	call dotvect(pn,d21,b(1))
	call dotvect(pm,d23,b(2))
	call dotvect(p2,u1,b(3))
c
	v1(1) = d21(1)
	v2(1) = d21(2)
	v3(1) = d21(3)
c
	v1(2) = d23(1)
	v2(2) = d23(2)
	v3(2) = d23(3)
c
	v1(3) = u1(1)
	v2(3) = u1(2)
	v3(3) = u1(3)
c
	call detvect(v1,v2,v3,det)
c
	call detvect(b,v2,v3,dd(1))
	call detvect(v1,b,v3,dd(2))
	call detvect(v1,v2,b,dd(3))
c
	do 100 i = 1,3
		v2(i) = dd(i)/det
100	continue
c
	call diffvect(p2,v2,v1)
	call normvect(u1,dist1)
	call normvect(v1,dist2)
c
	if(ip.eq.1) then
		alpha = -dsqrt(dist**2 - dist2**2)
	else
		alpha = dsqrt(dist**2 - dist2**2)
	endif
c
	do 200 i = 1,3
		u1(i) = (alpha/dist1) * u1(i)
200	continue
c
	call addvect(u1,v2,p4)
c
	return
	end
c	sp2.f		version 1 8/3/1990		Patrice Koehl
c
c	This subroutine adds a proton H4 on a carbon with a sp2
c	conformation
c
c		           C1		 
c			      \		 
c			        C2 ---- C3
c			      /	
c			   H4
c
c	Input of the program :
c				ang1	: angle H4-c2-c1
c				dist	: length of C2-H4 
c				p1	: C1
c				p2	: C2
c				p3	: C3
c	Output of the program :
c				p4	: H4 (or X)
c
	subroutine sp2(p1,p2,p3,p4,ang1,dist)
c
	real*8	p1(3),p2(3),p3(3),p4(3),ang1,ang2,dist
	real*8	v1(3),v2(3),v3(3),b(3),u1(3),dd(3),det
	real*8	dist1,dist2,d12(3),d23(3),dot,cos1,ang3,pi
c
	integer i
c
	pi = acos(-1.d0)
c
	call diffvect(p1,p2,d12)
	call diffvect(p2,p3,d23)
c
	call normvect(d12,dist1)
	call normvect(d23,dist2)
c
	call dotvect(d12,d23,dot)
	cos1 = -dot/(dist1*dist2)
	ang3 = acos(cos1)
	ang2 = 2*pi - ang1 - ang3
c
	b(1) = - dist1*dist*cos(ang1)
	b(2) = dist2*dist*cos(ang2)
	b(3) = 0.d0
c
	v1(1) = d12(1)
	v2(1) = d12(2)
	v3(1) = d12(3)
c
	v1(2) = d23(1)
	v2(2) = d23(2)
	v3(2) = d23(3)
c
	call crossvect(d12,d23,u1)
c
	v1(3) = u1(1)
	v2(3) = u1(2)
	v3(3) = u1(3)
c
	call detvect(v1,v2,v3,det)
c
	call detvect(b,v2,v3,dd(1))
	call detvect(v1,b,v3,dd(2))
	call detvect(v1,v2,b,dd(3))
c
	do 100 i = 1,3 
		u1(i) = dd(i)/det
100	continue
c
	call addvect(u1,p2,p4)
c
	return
	end
c	Buildfour.f		Version 1 8/2/1990	Patrice Koehl
c
c	This subroutine, knowing the position of 3 points and the torsion
c	angle that relates them to a 4th point, get the coordinates of
c	this fourth point
c	For details on how it works, see author
c	input of the program :
c			p1,p2,p3 and p4 are the four points considered;
c			tor is the torsional angle (in radian)
c			dist is the distance between p3 and p4
c			ang is ang(p2p3,p3p4)
c	output :
c			p4, coordinates of the fourth points
c
	subroutine buildfour(p1,p2,p3,p4,tor,ang,dist)
c
	integer i
	real*8	p1(3),p2(3),p3(3),p4(3),tor,dist,ang
	real*8	d21(3),d23(3),u1(3),u2(3)
	real*8	fact,b(3),v1(3),v2(3),v3(3),fact2
	real*8	norm,norm2,det,detx(3),sin1
	real*8	norm3,dot,cosm1
c
	call diffvect(p1,p2,d21)
	call diffvect(p2,p3,d23)
c
	call crossvect(d21,d23,u1)
	call dotvect(d21,d23,dot)
c
	call normvect(u1,norm)
	call normvect(d21,norm2)
	call normvect(d23,norm3)
c
	cosm1 = acos(dot/(norm2*norm3))
	sin1 = dabs(sin(cosm1))
c
	v1(1) = u1(1)/norm
	v2(1) = u1(2)/norm
	v3(1) = u1(3)/norm
c
	v1(2) = d21(1)/norm2
	v2(2) = d21(2)/norm2
	v3(2) = d21(3)/norm2
c
	v1(3) = d23(1)
	v2(3) = d23(2)
	v3(3) = d23(3)
c
	b(1)  = cos(tor)
	b(2)  = sin(tor)*sin1
	b(3)  = 0.d0
c
	call detvect(v1,v2,v3,det)
c
	call detvect(b,v2,v3,detx(1))
	call detvect(v1,b,v3,detx(2))
	call detvect(v1,v2,b,detx(3))
c
	do 100 i = 1,3
		u2(i) = detx(i)/det
100	continue
c
	fact = dist*dabs(sin(ang))/norm3
	fact2 = dist*dabs(cos(ang))/norm3
c
	call crossvect(u2,d23,v1)
c
	do 200 i = 1,3
		v1(i) = v1(i) * fact
		v2(i) = d23(i)* fact2
200	continue
c
	call addvect(v1,v2,v3)
	call addvect(v3,p3,p4)
c
	return
	end
c	Readpdb.f		Version 1 5/5/1992	Patrice Koehl
c
c	This program reads the present state of the molecule in
c	PDB format
c	At this stage, only one monomer of the molecule is considered
c
	subroutine readpdb(fname,keep_h,label)
c
	include 'pocket.h'
c
	real*8	x,y,z,ang,dist,pi,pif,torsion
	real*8  pa(3),pb(3),pc(3),pd(3),pe(3),pf(3)
	real*8	coord(ncortot),fulres(3*natrestot)
c
	integer	nseq,nat,i,natom,ierr,keep_h,nfulres,ip
	integer	nres,ntest,idum2,j
	integer	itype(nrestot),listatom(nrestot)
	integer	check(natot),checkat(natrestot)
	integer	chain(nrestot)
c
	character fname*30,atom(natot)*4,seq(nrestot)*4
	character label(natot)*54
	character name*4,name1*4,record*80
	character ch,ch1
c
	common /protein/ nseq,itype,chain,seq
	common /xyz/     coord,check,natom,listatom
	common /names/   atom
c
1	format(12x,a4,1x,4x,a1,i4,4x,3f8.3)
2	format(a)
c
	pi = acos(-1.d0)
	pif = pi/180.d0
c
	do 50 i = 1,natom
		check(i) = 0
50	continue
c
	open(unit=1,file=fname,status='old')
c
	ch1 = 'X'
	nres = 1
	ntest = -10
c
100     read(1,2,end=200) record
        if(record(1:6).ne.'ATOM  ') goto 100
	if(record(27:27).ne.' ') goto 100
        read(record,1) name1,ch,idum2,x,y,z
	if(idum2.ne.ntest) then
c		if(ntest.ne.-10.and.idum2-ntest.ne.1) then
c			write(6,*) 'Warning : non consecutive residues !!'
c			write(6,*) 'Residues : ',ntest,idum2
c		endif
		ntest = idum2
		if(nres.eq.1) ch1 = ch
		nres = nres + 1
		if(ch.ne.ch1) then
			ch1 = ch
			chain(nres) = 1
			chain(nres-1) = 2
		else
			chain(nres) = 0
		endif
	endif
        if(name1(1:1).eq.' ') then
                name(1:4) = name1(2:4)//' '
        else
                name = name1
        endif
	call atompos(nres,name,nat,ierr)
c	if(ierr.eq.0.or.ierr.eq.-1) goto 100
	if(check(nat).ne.0) goto 100
	check(nat) = check(nat) + 1
	coord(3*nat-2) = x
	coord(3*nat-1) = y
	coord(3*nat)   = z
	label(nat)     = record(1:54)
	goto 100
200	continue
c
	if(nres.ne.nseq-1) then
		write(6,*) 'inconsistence in residue number !!'
		stop
	endif
c
	chain(nseq) = 2
c
	close(unit=1)
c
	nat = 0
c
c	2. Check if there are missing N, CA, C and O
c
	nat = 2
	do 400 i = 2,nseq-1
		do 300 j = 1,listatom(i)
			nat = nat + 1
			if(j.eq.2.or.j.eq.4) goto 300
			if(check(nat).ne.1.and.j.le.6) then
				write(6,*) 'Problem in the PDB file !'
				write(6,*) 'Missing backbone atom : ',
     1				atom(nat),' in residue #',i,seq(i)
				stop
			endif
300		continue
400	continue
c
c	1. Build all missing CB atoms
c
	nat = listatom(1)
c
	do 1600 i = 2,nseq-1
c
		if(itype(i).eq.1.or.check(nat+4).eq.1) goto 1500
c
		do 1300 j = 1,3
			pa(j) = coord(3*nat+j)
			pb(j) = coord(3*(nat+2)+j)
			pc(j) = coord(3*(nat+4)+j)
1300		continue
c
		ip = 2
		dist = 1.53d0
		ang = 109.d0*pif
		call nsp3(pa,pb,pc,pd,ang,ang,ip,dist)
c
		do 1400 j = 1,3
			coord(3*(nat+3)+j) = pd(j)
1400		continue
c
1500		continue
c
		nat = nat + listatom(i)
c
1600	continue
c
c	2. Build Oxygen of Cter
c
	nat = nat - listatom(nseq-1)
c
	do 1700 i = 1,3
		pa(i) = coord(3*(nat+2)+i)
		pb(i) = coord(3*(nat+4)+i)
		pc(i) = coord(3*(nat+5)+i)
1700	continue
	dist = 1.20d0
	ang = 120.d0*pif
	call sp2(pa,pb,pc,pd,ang,dist)
c
	nat = nat + listatom(nseq-1)
	if(check(nat+1).eq.0) then
		do 1800 i = 1,3
			coord(3*nat+i) = pd(i)
1800		continue
	endif
c
	if(keep_h.eq.0) return
c
c	3. Build all missing backbone HN :
c
	nat = listatom(1)
c
c	2.a Build all three HN at the Nter of the protein, if missing
c
	do 500 i = 1,3
		pa(i) = coord(3*nat+i)
		pb(i) = coord(3*(nat+2)+i)
		pc(i) = coord(3*(nat+4)+i)
500	continue
	torsion	= pi
	ang 	= 111.d0*pif
	dist 	= 1.0d0
	call buildfour(pc,pb,pa,pd,torsion,ang,dist)
c
	ip = 2
	call nsp3(pb,pa,pd,pe,ang,ang,ip,dist)
	ip = 1
	call nsp3(pb,pa,pd,pf,ang,ang,ip,dist)
c
	if(check(1).eq.0.or.check(2).eq.0.or.check(4).eq.0) then
		do 600 i = 1,3
			coord(i) 	= pd(i)
			coord(i+3) 	= pe(i)
			coord(i+9) 	= pf(i)
600		continue
	endif
c
	nat = nat + listatom(2)
	do 700 i = 1,3
		pa(i) = pc(i)
700	continue
c
c	2.b for all residues, check if HN was read in the PDB file,
c	otherwise, build it through standard geometry
c
	do 1200 i = 3,nseq-1
c
		if(check(nat+2).ne.0) then
			do 750 j = 1,3
				pb(j) = coord(3*nat+j)
				pd(j) = coord(3*(nat+1)+j)
750			continue
			goto 1000
		endif
c
		do 800 j = 1,3
			pb(j) = coord(3*nat+j)
			pc(j) = coord(3*(nat+2)+j)
800		continue
c
		ang = 123.d0*pif
		call sp2(pa,pb,pc,pd,ang,dist)
c
		do 900 j = 1,3
			coord(3*(nat+1)+j) = pd(j)
900		continue
c
1000		continue
c
		do 1100 j = 1,3
			pa(j) = coord(3*(nat+4)+j)
1100		continue
c
		nat = nat + listatom(i)
c
1200	continue
c
c	4. Build missing polar H on sidechains
c
	nat = listatom(1)
c
	do 2200 i = 2,nseq-1
c
		nfulres = listatom(i)
		do 1900 j = 1,3*nfulres
			fulres(j) = coord(3*nat + j)
1900		continue
		do 2000 j = 1,nfulres
			checkat(j) = check(nat+j)
2000		continue
c
		call addh(itype(i),fulres,nfulres,checkat)
c
		do 2100 j = 1,3*nfulres
			coord(3*nat+j) = fulres(j)
2100		continue
c
		nat = nat + listatom(i)
c
2200	continue
c
c
	return
	end
c	Vector
c
c	This file contains several subroutines that can be used for any
c	vector operations (vector in 3D cartesian space)
c
c	This includes :	crossvect	: cross vector of two vectors
c			dotvect		: dot product of two vectors
c			normvect	: norm of a vector
c			detvect		: determinant of three vectors
c			diffvect	: substract two vectors
c			addvect		: add two vectors
c
c	For each subroutine : u is for vectors (arrays of size 3)
c			      all other value are scalar
c			      calculations are done in double precision
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
c	1 . crossvect :
c
	subroutine crossvect(u1,u2,u3)
c
	real*8	u1(3),u2(3),u3(3)
c
	u3(1) = u1(2)*u2(3) - u1(3)*u2(2)
	u3(2) = -u1(1)*u2(3) + u1(3)*u2(1)
	u3(3) = u1(1)*u2(2) - u1(2)*u2(1)
c
	return
	end
c
c	2. dotvect :
c
	subroutine dotvect(u1,u2,dot)
c
	integer	i
c
	real*8	u1(3),u2(3),dot
c
	dot = 0.d0
	do 100 i = 1,3
		dot = dot + u1(i)*u2(i)
100	continue
c
	return
	end
c
c	3. normvect :
c
	subroutine normvect(u1,norm)
c
	real*8	u1(3),norm
c
	call dotvect(u1,u1,norm)
	norm = dsqrt(norm)
c
	return
	end
c
c	4. detvect :
c
	subroutine detvect(u1,u2,u3,det)
c
	real*8	u1(3),u2(3),u3(3),det,u4(3)
c
	call crossvect(u2,u3,u4)
	call dotvect(u1,u4,det)
c
	return
	end
c
c	5. diffvect :
c
	subroutine diffvect(u1,u2,u3)
c
	real*8	u1(3),u2(3),u3(3)
c
	integer i
c
	do 100 i = 1,3
		u3(i) = u2(i) - u1(i)
100	continue
c
	return
	end
c
c	6. addvect :
c
	subroutine addvect(u1,u2,u3)
c
	real*8	u1(3),u2(3),u3(3)
c
	integer i
c
	do 100 i = 1,3
		u3(i) = u1(i) + u2(i)
100	continue
c
	return
	end
c
c	7. Normalise a vector : given a vector u1, output u1/norm(u1) :
c
	subroutine unitvector(u1,u2)
c
	real*8  u1(3),u2(3),norm
c
	integer i
c
	call normvect(u1,norm)
c
	do 100 i = 1,3
		u2(i) = u1(i)/norm
100	continue
c
	return
	end
