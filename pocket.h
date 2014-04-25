C	Pocket.h	Version 1 1/16/2002	Patrice Koehl
c
c	This include file contains all the dimensions of the
c	arrays needed for computing the triangulation
c
	integer	nresdef,nrestot,nangtot,natot,ncortot
	integer	natom_type,natrestot,naspmax
	integer	npointmax,ntetra_max
	integer	nfacet_max,nfreemax
	integer	ntrig_max,nedge_max
	integer	nvoid_max,nmouth_max
	integer nlistmax,ntet_listmax
c
c	Definitions:
c
c		nresdef		: maximum number of types of residues
c		nrestot		: maximum number of residue in protein
c		nangtot		: maximum number of dihedral angle in 1 protein
c		natot		: maximum number of atom in 1 protein
c		ncortot		: 3 * natot
c		natom_type	: maximum number of atom types
c		natrestot	: maximum number of atoms in one residue
c		naspmax		: maximum number of types of atoms for solvation
c		npointmax	: maximum number of points for triangulation
c		ntetra_max 	: maximum number of tetrahedra in DT
c		ntrig_max	: maximum number of triangles in DT
c		nedge_max	: maximum number of edges in DT
c		nfacet_max	: maximum number of facets to be flipped
c				  (used for DT only)
c		nfreemax	: maximum number of "free" tetrahedra after
c				  each flipping sequence. These space are
c				  then re-used.
c
	parameter	(nresdef	= 25)
	parameter	(nrestot	= 5000)
	parameter	(nangtot	= 40000)
	parameter	(natot		= 100000)
	parameter	(ncortot	= 3*natot)
	parameter	(natom_type	= 30)
	parameter	(natrestot	= 20)
	parameter	(naspmax	= 8)
	parameter	(npointmax	= 100000)
	parameter	(ntetra_max	= 500000)
	parameter	(ntrig_max	= 1000000)
	parameter	(nedge_max	= 700000)
	parameter	(nfacet_max	= 5000)
	parameter	(nfreemax	= 3000)
	parameter	(nvoid_max	= 10000)
	parameter	(nlistmax	= 100000)
	parameter	(ntet_listmax   = 500)
	parameter	(nmouth_max	= 10000)
c
