/*	tetra_flow_gmp.c	Version 1 11/27/2000	Patrice Koehl     */
/*									  */
/*  This is the C version of tetra_flow.f, which performs all operations  */
/*  with multi precision arithmetics, using the package GMP               */
/*									  */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* Copyright (C) 2002 Patrice Koehl                                       */
/*                                                                        */
/* This library is free software; you can redistribute it and/or          */
/* modify it under the terms of the GNU Lesser General Public             */
/* License as published by the Free Software Foundation; either           */
/* version 2.1 of the License, or (at your option) any later version.     */
/*                                                                        */
/* This library is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of         */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU       */
/* Lesser General Public License for more details.                        */
/*                                                                        */
/* You should have received a copy of the GNU Lesser General Public       */
/* License along with this library; if not, write to the Free Software    */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA*/
/*                                                                        */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* Includes :								  */

#include <stdio.h>
#include <math.h>
#include "gmp.h"

/*------------------------------------------------------------------------*/
/*Local procedures:							  */

void transfer_to_mp(double *a, double *b, double *c, double *d,
		double *ra, double *rb, double *rc, double *rd,
		double *scale, mpz_t *a_mp, mpz_t wa, mpz_t *b_mp,
		mpz_t wb, mpz_t *c_mp, mpz_t wc, mpz_t *d_mp, mpz_t wd);

void check_facet_gmp(mpz_t *S, mpz_t De1,
		mpz_t De2, mpz_t De3, mpz_t Dmnpq, mpz_t Dmnp,
		int *testa);

void tetra_flow_gmp_(double *a, double *b, double *c, double *d,
		double *ra, double *rb, double *rc, double *rd,
		int *test, double *scale);

/*------------------------------------------------------------------------*/
/* transfer_to_mp:                                                  
	This procedure transfers the coordinates of the 4 vertices of
	a tetrahedron in multi precision format                           */

void transfer_to_mp(double *a, double *b, double *c, double *d,
		double *ra, double *rb, double *rc, double *rd,
		double *scale, mpz_t *a_mp, mpz_t wa, mpz_t *b_mp,
		mpz_t wb, mpz_t *c_mp, mpz_t wc, mpz_t *d_mp, mpz_t wd)

{

	int i,ivalue;
	double value;
	mpz_t temp1,temp2;

	mpz_init(temp1); mpz_init(temp2);

	for ( i=0; i<3 ; i++)
	{
		value = a[i]*(*scale); ivalue = (int) floor(value); 
		mpz_set_si(a_mp[i+1],ivalue);
		value = b[i]*(*scale); ivalue = (int) floor(value); 
		mpz_set_si(b_mp[i+1],ivalue);
		value = c[i]*(*scale); ivalue = (int) floor(value); 
		mpz_set_si(c_mp[i+1],ivalue);
		value = d[i]*(*scale); ivalue = (int) floor(value); 
		mpz_set_si(d_mp[i+1],ivalue);
	}
	value=(*ra)*(*scale);ivalue=(int) floor(value); 
	mpz_set_si(temp1,ivalue); mpz_mul(temp1,temp1,temp1);
	mpz_mul(temp2,a_mp[3],a_mp[3]), mpz_sub(temp1,temp2,temp1);
	mpz_mul(temp2,a_mp[2],a_mp[2]), mpz_add(temp1,temp2,temp1);
	mpz_mul(temp2,a_mp[1],a_mp[1]), mpz_add(wa,temp2,temp1);
	value=(*rb)*(*scale);ivalue=(int) floor(value); 
	mpz_set_si(temp1,ivalue); mpz_mul(temp1,temp1,temp1);
	mpz_mul(temp2,b_mp[3],b_mp[3]), mpz_sub(temp1,temp2,temp1);
	mpz_mul(temp2,b_mp[2],b_mp[2]), mpz_add(temp1,temp2,temp1);
	mpz_mul(temp2,b_mp[1],b_mp[1]), mpz_add(wb,temp2,temp1);
	value=(*rc)*(*scale);ivalue=(int) floor(value); 
	mpz_set_si(temp1,ivalue); mpz_mul(temp1,temp1,temp1);
	mpz_mul(temp2,c_mp[3],c_mp[3]), mpz_sub(temp1,temp2,temp1);
	mpz_mul(temp2,c_mp[2],c_mp[2]), mpz_add(temp1,temp2,temp1);
	mpz_mul(temp2,c_mp[1],c_mp[1]), mpz_add(wc,temp2,temp1);
	value=(*rd)*(*scale);ivalue=(int) floor(value);
	mpz_set_si(temp1,ivalue); mpz_mul(temp1,temp1,temp1);
	mpz_mul(temp2,d_mp[3],d_mp[3]), mpz_sub(temp1,temp2,temp1);
	mpz_mul(temp2,d_mp[2],d_mp[2]), mpz_add(temp1,temp2,temp1);
	mpz_mul(temp2,d_mp[1],d_mp[1]), mpz_add(wd,temp2,temp1);

	mpz_clear(temp1); mpz_clear(temp2);
}

/*------------------------------------------------------------------------*/
/* check_facet_gmp:
   This program checks if a face of a tetrahedron is attached to the fourth
   vertex of the tetrahedron

	Input:

	For the three points m,n,p that form the triangles, the program
	needs as input the following determinants:

	S(i+j-2) = Minor(m,n,p,i,j,0)= det | m(i)  m(j)  1 |
					   | n(i)  n(j)  1 |
					   | p(i)  p(j)  1 |

If q is the fourth point of the tetrahedron,

	De1 = Minor(m,n,p,q,2,3,4,0)
	De2 = Minor(m,n,p,q,1,3,4,0)
	De3 = Minor(m,n,p,q,1,2,4,0)
	Dmnpq= Minor(m,n,p,q,1,2,3,0)
	Dmnp  = Minor(m,n,p,1,2,3)

	Output:

	testa	: flag set to 1 if the fourth point q of the tetrahedron
		  that contains the triangle (m,n,p) in inside the
		  circumsphere of (m,n,p)
		  This is equivalent to testing if the orthocenter C of (m,n,p,q)
		  is on the same side of (m,n,p) than q or not
		  (if q is inside the circumsphere of (m,n,p), C and q are
		  on opposite side)

*/

void check_facet_gmp(mpz_t *S, mpz_t De1, 
	mpz_t De2, mpz_t De3, mpz_t Dmnpq, mpz_t Dmnp,int *testa)

{
	int i,coef,icomp;


/*	Local GMP variables */

	mpz_t temp1, temp2, temp3;
	mpz_t dtest;

/* 	Initialise local GMP variables */

	mpz_init(temp1); mpz_init(temp2); mpz_init(temp3);
	mpz_init (dtest);

/* First check if triangle is attached                           */

	mpz_set_si(temp1,0);
	for (i=1; i<4; i++)
	{
		mpz_mul(temp2,S[i],S[i]);
		mpz_add(temp1,temp1,temp2);
	}

	if(!(*testa))
	{
		mpz_mul(temp2,Dmnpq,Dmnp);
		coef = -2;
		mpz_mul_si(temp2,temp2,coef);
		mpz_mul(temp3,De3,S[1]);
		mpz_add(temp2,temp3,temp2);
		mpz_mul(temp3,De2,S[2]);
		mpz_add(temp2,temp3,temp2);
		mpz_mul(temp3,De1,S[3]);
		mpz_add(temp2,temp3,temp2);
		mpz_mul(dtest,temp1,temp2);

		icomp = mpz_sgn(dtest);
		if(icomp > 0) {
			(*testa)=1;
		}
		else if (icomp == 0) {
			printf("Warning: degeneracy in flow_gmp\n");
			(*testa)=0;
		}
		else {
			(*testa)=0;
		}
	}

/* 	Clear local GMP variables */

	mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3);
	mpz_clear (dtest);
}
/*------------------------------------------------------------------------*/
/* tetra_flow_gmp_	Version 1 11/24/2000	Patrice Koehl
 
	This subroutine computes the flow from a tetrahedron [A,B,C,D]

	The "flow" at a face (I,J,K) is defined based on the position
	of the orthocenter CO of the tetrahedron with respect to the plane
	(I,J,K): if CO and L (where L is the fourth point of the tetrahedron)
	are on the same side the flow is "IN", otherwise the flow is "OUT"
	
 	This position test is in fact equivalent to the "attached" test used
	for building the Alpha complex: the "attached" test checks if any 
	fourth point L (in A, B, C, D) of the tetrahedron is "hidden" by its 
	opposite face [I,J,K], i.e. is interior to the cicumsphere of [I,J,K]
	(if L is inside the circumsphere of (I,J,K), CO and L are
		  on opposite side)
 
 	Computation extends to all four faces of the tetrehedron
 
	This procedure works with Multiple Precision Integer Arithmetics  (MPIA)
 	The package GMP is used for MPIA (with a C wrapper)

*/ 
void tetra_flow_gmp_(double *a, double *b, double *c, double *d,
		double *ra, double *rb, double *rc, double *rd,
		int *test, double *scale)

{
	int i,j,k,coef,test_a;
	int ivalue;
	double value;

/*	Local GMP variables	*/

	mpz_t wa,wb,wc,wd;
	mpz_t temp1,temp2,temp3;
	mpz_t val1,val2,val3,val4;
	mpz_t Dabc,Dabd,Dacd,Dbcd;
	mpz_t Det1,Det2,Det3,Det4;

	mpz_t a_mp[4], b_mp[4],c_mp[4],d_mp[4];
	mpz_t Sab[4],Sac[4],Sad[4],Sbc[4],Sbd[4],Scd[4];
	mpz_t Sa[4],Sb[4],Sc[4],Sd[4];
	mpz_t Deter[4];

/*	Initialise local GMP variables */

	for (i = 0; i < 4; i++) 
	{
		mpz_init(a_mp[i]);
		mpz_init(b_mp[i]);
		mpz_init(c_mp[i]);
		mpz_init(d_mp[i]);
	}
	mpz_init(wa);mpz_init(wb);mpz_init(wc);mpz_init(wd);

	mpz_init(temp1); mpz_init(temp2); mpz_init(temp3);
	mpz_init(val1); mpz_init(val2); mpz_init(val3); mpz_init(val4);

	for (i=0; i < 4; i++)
	{
		mpz_init(Sab[i]); mpz_init(Sac[i]); mpz_init(Sad[i]);
		mpz_init(Sbc[i]); mpz_init(Sbd[i]); mpz_init(Scd[i]);
		mpz_init(Sa[i]); mpz_init(Sb[i]); mpz_init(Sc[i]); 
		mpz_init(Sd[i]);
		mpz_init(Deter[i]);
	}

	mpz_init(Dabc); mpz_init(Dabd); mpz_init(Dacd); mpz_init(Dbcd);
	mpz_init(Det1); mpz_init(Det2); mpz_init(Det3); 
	mpz_init(Det4);

/*	Transfer data in multiple precision */

	transfer_to_mp(a,b,c,d,ra,rb,rc,rd,scale,a_mp,wa,
	b_mp,wb,c_mp,wc,d_mp,wd);

/*	1. Computes all Minors Smn(i+j-2)= M(m,n,i,j) = Det | m(i)  m(j) |
						            | n(i)  n(j) |
	for all i in [1,2] and all j in [i+1,3]                         */

	for (i=1;  i<3; i++)
	{
		for (j=i+1; j<4 ; j++)
		{
			k=i+j-2;
			mpz_mul(temp1,a_mp[j],b_mp[i]); 
			mpz_mul(temp2,a_mp[i],b_mp[j]);
			mpz_sub(Sab[k],temp2,temp1);
			mpz_mul(temp1,a_mp[j],c_mp[i]); 
			mpz_mul(temp2,a_mp[i],c_mp[j]);
			mpz_sub(Sac[k],temp2,temp1);
			mpz_mul(temp1,a_mp[j],d_mp[i]); 
			mpz_mul(temp2,a_mp[i],d_mp[j]);
			mpz_sub(Sad[k],temp2,temp1);
			mpz_mul(temp1,b_mp[j],c_mp[i]); 
			mpz_mul(temp2,b_mp[i],c_mp[j]);
			mpz_sub(Sbc[k],temp2,temp1);
			mpz_mul(temp1,b_mp[j],d_mp[i]); 
			mpz_mul(temp2,b_mp[i],d_mp[j]);
			mpz_sub(Sbd[k],temp2,temp1);
			mpz_mul(temp1,c_mp[j],d_mp[i]); 
			mpz_mul(temp2,c_mp[i],d_mp[j]);
			mpz_sub(Scd[k],temp2,temp1);
		}
	}

/*	Now compute all Minors 
		Sq(i+j-2) = M(m,n,p,i,j,0) = Det | m(i) m(j) 1 |
		       			         | n(i) n(j) 1 |
						 | p(i) p(j) 1 |

	and all Minors
		Det(i+j-2) = M(m,n,p,q,i,j,4,0) = Det | m(i) m(j) m(4) 1 |
						      | n(i) n(j) n(4) 1 |
						      | p(i) p(j) p(4) 1 |
						      | q(i) q(j) q(4) 1 |

	m,n,p,q are the four vertices of the tetrahedron, i and j correspond
	to two of the coordinates of the vertices, and m(4) refers to the
	"weight" of vertices m                                           */
 
	for (i=1; i<4; i++)
	{
		mpz_sub(temp1,Scd[i],Sbd[i]); mpz_add(Sa[i],temp1,Sbc[i]);
		mpz_mul(temp2,Sa[i],wa);
		mpz_sub(temp1,Scd[i],Sad[i]); mpz_add(Sb[i],temp1,Sac[i]);
		mpz_mul(temp3,Sb[i],wb); mpz_sub(temp2,temp2,temp3);
		mpz_sub(temp1,Sbd[i],Sad[i]); mpz_add(Sc[i],temp1,Sab[i]);
		mpz_mul(temp3,Sc[i],wc); mpz_add(temp2,temp2,temp3);
		mpz_sub(temp1,Sbc[i],Sac[i]); mpz_add(Sd[i],temp1,Sab[i]);
		mpz_mul(temp3,Sd[i],wd); mpz_sub(Deter[i],temp2,temp3);
	}
 
/*
	Now compute the determinant needed to compute the radius of the
	circumsphere of the tetrahedron :

		Det1 = Minor(a,b,c,d,4,2,3,0)
		Det2 = Minor(a,b,c,d,1,3,4,0)
		Det3 = Minor(a,b,c,d,1,2,4,0)
		Det4 = Minor(a,b,c,d,1,2,3,0)
									*/

	mpz_set(Det1,Deter[3]);
	mpz_set(Det2,Deter[2]);
	mpz_set(Det3,Deter[1]);

	mpz_mul(temp1,a_mp[1],Sa[3]);mpz_mul(temp2,b_mp[1],Sb[3]);
	mpz_sub(temp3,temp1,temp2);
	mpz_mul(temp1,c_mp[1],Sc[3]);mpz_mul(temp2,d_mp[1],Sd[3]);
	mpz_sub(temp1,temp1,temp2);
	mpz_add(Det4,temp1,temp3);

/*
	Now compute all minors:
		Dmnp = Minor(m,n,p,1,2,3) = Det | m(1) m(2) m(3) |
						| n(1) n(2) n(3) |
						| p(1) p(2) p(3) |
									*/

	mpz_mul(temp1,a_mp[1],Sbc[3]); mpz_mul(temp2,b_mp[1],Sac[3]);
	mpz_sub(temp3,temp1,temp2);
	mpz_mul(temp1,c_mp[1],Sab[3]);mpz_add(Dabc,temp3,temp1);

	mpz_mul(temp1,a_mp[1],Sbd[3]); mpz_mul(temp2,b_mp[1],Sad[3]);
	mpz_sub(temp3,temp1,temp2);
	mpz_mul(temp1,d_mp[1],Sab[3]);mpz_add(Dabd,temp3,temp1);

	mpz_mul(temp1,a_mp[1],Scd[3]); mpz_mul(temp2,c_mp[1],Sad[3]);
	mpz_sub(temp3,temp1,temp2);
	mpz_mul(temp1,d_mp[1],Sac[3]);mpz_add(Dacd,temp3,temp1);

	mpz_mul(temp1,b_mp[1],Scd[3]); mpz_mul(temp2,c_mp[1],Sbd[3]);
	mpz_sub(temp3,temp1,temp2);
	mpz_mul(temp1,d_mp[1],Sbc[3]);mpz_add(Dbcd,temp3,temp1);


/* 
 	Now check all four faces of the tetrahedra:
 	We look if the fourth vertex is "hidden" by the face of 
 	interest
 
 									*/
	mpz_neg(val1,Det1); mpz_neg(val2,Det2);
	mpz_neg(val3,Det3); mpz_neg(val4,Det4);

/* First check face abc */

	if(!test[3])
	{
		test_a = test[3];
		check_facet_gmp(Sd,Det1,Det2,Det3,Det4,Dabc,
			&test_a);
		test[3]=test_a;
	}

/* Now check face abd */
 
	if(!test[2])
	{
		test_a = test[2];
		check_facet_gmp(Sc,val1,val2,val3,val4,Dabd,
			&test_a);
		test[2]=test_a;
	}

/* Now check face acd */
 
	if(!test[1])
	{
		test_a = test[1];
		check_facet_gmp(Sb,Det1,Det2,Det3,Det4,Dacd,
			&test_a);
		test[1]=test_a;
	}
 
/* Now check face bcd */
 
	if(!test[0])
	{
		test_a = test[0];
		check_facet_gmp(Sa,val1,val2,val3,val4,Dbcd,
			&test_a);
		test[0]=test_a;
	}

/*	Clear local GMP variables */

	for (i = 0; i < 4; i++) 
	{
		mpz_clear(a_mp[i]);
		mpz_clear(b_mp[i]);
		mpz_clear(c_mp[i]);
		mpz_clear(d_mp[i]);
	}
	mpz_clear(wa);mpz_clear(wb);mpz_clear(wc);mpz_clear(wd);
	mpz_clear(val1); mpz_clear(val2); mpz_clear(val3); mpz_clear(val4);

	mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3);

	for (i=0; i < 4; i++)
	{
		mpz_clear(Sab[i]); mpz_clear(Sac[i]); mpz_clear(Sad[i]);
		mpz_clear(Sbc[i]); mpz_clear(Sbd[i]); mpz_clear(Scd[i]);
		mpz_clear(Sa[i]); mpz_clear(Sb[i]); mpz_clear(Sc[i]); 
		mpz_clear(Sd[i]);
		mpz_clear(Deter[i]);
	}

	mpz_clear(Dabc); mpz_clear(Dabd); mpz_clear(Dacd); mpz_clear(Dbcd);
	mpz_clear(Det1); mpz_clear(Det2); mpz_clear(Det3);
	mpz_clear(Det4);

}
