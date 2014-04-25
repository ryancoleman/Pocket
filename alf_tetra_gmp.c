/*	alf_tetra_gmp.c	Version 1 11/27/2000	Patrice Koehl             */
/*									  */
/*  This is the C version of alf_tetra.f, which performs all operations   */
/*  with multi precision arithmetics, using the package GMP               */
/*									  */
/*------------------------------------------------------------------------*/
/*									  */
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
/*									  */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* Includes :								  */

#include <stdio.h>
#include <math.h>
#include "gmp.h"

/*------------------------------------------------------------------------*/
/*Local procedures:							  */

void transfer_coord_to_mp(double *a, double *b, double *c, double *d,
		double *ra, double *rb, double *rc, double *rd,
		double *scale, mpz_t *a_mp, mpz_t wa, mpz_t *b_mp,
		mpz_t wb, mpz_t *c_mp, mpz_t wc, mpz_t *d_mp, mpz_t wd);

void check_edge_gmp(mpz_t *m_mp, mpz_t *n_np, mpz_t wm, mpz_t wn,
		mpz_t *Dab, mpz_t *Sab, mpz_t *Sc, mpz_t *Sd, 
		mpz_t *Tc, mpz_t *Td, int *testr, int *testa, mpz_t alp);

void check_triangle_gmp(mpz_t *S, mpz_t *T, mpz_t *U, mpz_t De1,
		mpz_t De2, mpz_t De3, mpz_t Dmnpq, mpz_t Dmnp,
		int *testr, int *testa, mpz_t alp);

void alf_tetra_gmp_(double *a, double *b, double *c, double *d,
		double *ra, double *rb, double *rc, double *rd,
		int *testr, int *testa, double *scale, double *alpha);

/*------------------------------------------------------------------------*/
/* transfer_coord_to_mp:                                                  
	This procedure transfers the coordinates of the 4 vertices of
	a tetrahedron in multi precision format                           */

void transfer_coord_to_mp(double *a, double *b, double *c, double *d,
		double *ra, double *rb, double *rc, double *rd,
		double *scale, mpz_t *a_mp, mpz_t wa, mpz_t *b_mp,
		mpz_t wb, mpz_t *c_mp, mpz_t wc, mpz_t *d_mp, mpz_t wd)

{

	int i;
	double value;
	mpz_t temp1,temp2;

	mpz_init(temp1); mpz_init(temp2);

	for ( i=0; i<3 ; i++)
	{
		value = a[i]*(*scale); 
		mpz_set_d(a_mp[i+1],value);
		value = b[i]*(*scale);  
		mpz_set_d(b_mp[i+1],value);
		value = c[i]*(*scale); 
		mpz_set_d(c_mp[i+1],value);
		value = d[i]*(*scale); 
		mpz_set_d(d_mp[i+1],value);
	}
	value=(*ra)*(*scale);
	mpz_set_d(temp1,value); mpz_mul(temp1,temp1,temp1);
	mpz_mul(temp2,a_mp[3],a_mp[3]), mpz_sub(temp1,temp2,temp1);
	mpz_mul(temp2,a_mp[2],a_mp[2]), mpz_add(temp1,temp2,temp1);
	mpz_mul(temp2,a_mp[1],a_mp[1]), mpz_add(wa,temp2,temp1);
	value=(*rb)*(*scale);
	mpz_set_d(temp1,value); mpz_mul(temp1,temp1,temp1);
	mpz_mul(temp2,b_mp[3],b_mp[3]), mpz_sub(temp1,temp2,temp1);
	mpz_mul(temp2,b_mp[2],b_mp[2]), mpz_add(temp1,temp2,temp1);
	mpz_mul(temp2,b_mp[1],b_mp[1]), mpz_add(wb,temp2,temp1);
	value=(*rc)*(*scale);
	mpz_set_d(temp1,value); mpz_mul(temp1,temp1,temp1);
	mpz_mul(temp2,c_mp[3],c_mp[3]), mpz_sub(temp1,temp2,temp1);
	mpz_mul(temp2,c_mp[2],c_mp[2]), mpz_add(temp1,temp2,temp1);
	mpz_mul(temp2,c_mp[1],c_mp[1]), mpz_add(wc,temp2,temp1);
	value=(*rd)*(*scale);
	mpz_set_d(temp1,value); mpz_mul(temp1,temp1,temp1);
	mpz_mul(temp2,d_mp[3],d_mp[3]), mpz_sub(temp1,temp2,temp1);
	mpz_mul(temp2,d_mp[2],d_mp[2]), mpz_add(temp1,temp2,temp1);
	mpz_mul(temp2,d_mp[1],d_mp[1]), mpz_add(wd,temp2,temp1);

	mpz_clear(temp1); mpz_clear(temp2);
}

/*------------------------------------------------------------------------*/
/* check_edge_gmp:
	This subroutine check if an edge of a tetrahedron is "attached"
	to one of the two remaining vertices of the tetrahedron (i.e.
	if one of the two vertices belong to the smallest circumsphere
	of the edge). For that, it needs:

	Input:
		a,b	: coordinate of the two vertices defining the edge
		wa,wb	: weights of these two vertices
		Dab	: minor(a,b,i,0) for all i=1,2,3
		Sab	: minor(a,b,i,j) for i = 1,2 and j =i+1,3
		Sc	: minor(a,b,c,i,j,0) for i=1,2 and j = i+1,3
		Sd	: minor(a,b,d,i,j,0) for i=1,2 and j = i+1,3
			  c and d are the two other vertices of the
			  tetrahedron
		Tc	: minor(a,b,c,i,4,0) for i = 1,2,3
		Td	: minor(a,b,d,i,4,0) for i = 1,2,3
		alpha	: value of alpha considered
		eps	: precision: if a floating point test lead to a
			  value below this precision, computation
			  switches to LIA
	Ouput:
		testr	: flag that defines if edge belongs to the
			  Alpha complex
		testa	: flag that defines if edge is attached or not
									  
For comments, see alpha_fixed.f which contains the fortran equivalence of this
routine
									*/

void check_edge_gmp(mpz_t *m_mp, mpz_t *n_mp, mpz_t wm, mpz_t wn,
		mpz_t *Dab, mpz_t *Sab, mpz_t *Sc, mpz_t *Sd, 
		mpz_t *Tc, mpz_t *Td, int *testr, int *testa, mpz_t alp)

{
	int i,j,coef;

/*	Local GMP variables */

	mpz_t temp1,temp2, temp3;
	mpz_t r_11, r_22, r_33, r_14, r_313, r_212, diff;
	mpz_t d0,d1,d2,d3,d4;
	mpz_t num,den,dtest;
	mpz_t res[4][5],res2_c[4][5],res2_d[4][5];

/* 	Initialise local GMP variables */

	mpz_init(temp1); mpz_init(temp2); mpz_init(temp3);
	mpz_init (r_11); mpz_init (r_22); mpz_init (r_33); mpz_init (r_14);
	mpz_init (r_313); mpz_init (r_212); mpz_init (diff);
	mpz_init (d0); mpz_init (d1); mpz_init (d2); mpz_init (d3);
	mpz_init (d4); 
	mpz_init (num); mpz_init (den); mpz_init (dtest);

	for (i= 0; i < 4; i++)
	{
		for (j=0; j < 5; j++)
		{
			mpz_init(res[i][j]);
			mpz_init(res2_c[i][j]);
			mpz_init(res2_d[i][j]);
		}
	}


/*	This is the "hidden1" part */

	mpz_sub(res[0][4],wm,wn);


	if( mpz_cmp(m_mp[1],n_mp[1]) != 0) 
	{
		for (i = 1; i < 4 ; i++)
		{
			mpz_set(res[0][i],Dab[i]);
			mpz_set(res2_c[i][4],Tc[i]);
			mpz_set(res2_d[i][4],Td[i]);
			mpz_mul(temp1,n_mp[i],wm);
			mpz_mul(temp2,m_mp[i],wn);
			mpz_sub(res[i][4],temp2,temp1);
		}
		mpz_set(res[1][2],Sab[1]); mpz_set(res[1][3],Sab[2]);
		mpz_set(res[2][3],Sab[3]);
		mpz_set(res2_c[1][2],Sc[1]); mpz_set(res2_c[1][3],Sc[2]);
		mpz_set(res2_c[2][3],Sc[3]);
		mpz_set(res2_d[1][2],Sd[1]); mpz_set(res2_d[1][3],Sd[2]);
		mpz_set(res2_d[2][3],Sd[3]);
	}
	else if ( mpz_cmp(m_mp[2],n_mp[2]) != 0)
	{
		mpz_set(res[0][1],Dab[2]); mpz_set(res[0][2],Dab[3]);
		mpz_set(res[0][3],Dab[1]);
		mpz_set(res[1][2],Sab[3]);
		mpz_neg(res[1][3],Sab[1]); mpz_neg(res[2][3],Sab[2]);
		mpz_mul(temp1,m_mp[2],wn); mpz_mul(temp2,n_mp[2],wm);
		mpz_sub(res[1][4],temp1,temp2);
		mpz_mul(temp1,m_mp[3],wn); mpz_mul(temp2,n_mp[3],wm);
		mpz_sub(res[2][4],temp1,temp2);
		mpz_mul(temp1,m_mp[1],wn); mpz_mul(temp1,n_mp[1],wm);
		mpz_sub(res[3][4],temp1,temp2);
		mpz_set(res2_c[1][2],Sc[3]);
		mpz_neg(res2_c[1][3],Sc[1]); mpz_neg(res2_c[2][3],Sc[2]);
		mpz_set(res2_d[1][2],Sd[3]);
		mpz_neg(res2_d[1][3],Sd[1]); mpz_neg(res2_d[2][3],Sd[2]);
		mpz_set(res2_c[1][4],Tc[2]); mpz_set(res2_c[2][4],Tc[3]);
		mpz_set(res2_c[3][4],Tc[1]);
		mpz_set(res2_d[1][4],Td[2]); mpz_set(res2_d[2][4],Td[3]);
		mpz_set(res2_d[3][4],Td[1]);
	}
	else if (  mpz_cmp(m_mp[3],n_mp[3]) != 0)
	{
		mpz_set(res[0][1],Dab[3]); mpz_set(res[0][2],Dab[1]);
		mpz_set(res[0][3],Dab[2]);
		mpz_neg(res[1][2],Sab[2]);
		mpz_neg(res[1][3],Sab[3]); mpz_set(res[2][3],Sab[1]);
		mpz_mul(temp1,m_mp[3],wn); mpz_mul(temp2,n_mp[3],wm);
		mpz_sub(res[1][4],temp1,temp2);
		mpz_mul(temp1,m_mp[1],wn); mpz_mul(temp2,n_mp[1],wm);
		mpz_sub(res[2][4],temp1,temp2);
		mpz_mul(temp1,m_mp[2],wn); mpz_mul(temp1,n_mp[2],wm);
		mpz_sub(res[3][4],temp1,temp2);
		mpz_neg(res2_c[1][2],Sc[2]);
		mpz_neg(res2_c[1][3],Sc[3]); mpz_set(res2_c[2][3],Sc[1]);
		mpz_neg(res2_d[1][2],Sd[2]);
		mpz_neg(res2_d[1][3],Sd[3]); mpz_set(res2_d[2][3],Sd[1]);
		mpz_set(res2_c[1][4],Tc[3]); mpz_set(res2_c[2][4],Tc[1]);
		mpz_set(res2_c[3][4],Tc[2]);
		mpz_set(res2_d[1][4],Td[3]); mpz_set(res2_d[2][4],Td[1]);
		mpz_set(res2_d[3][4],Td[2]);
	}
	else
	{
		mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3);
		mpz_clear (r_11); mpz_clear (r_22); mpz_clear (r_33); mpz_clear (r_14);
		mpz_clear (r_313); mpz_clear (r_212); mpz_clear (diff);
		mpz_clear (d0); mpz_clear (d1); mpz_clear (d2); mpz_clear (d3);
		mpz_clear (d4); 
		mpz_clear (num); mpz_clear (den); mpz_clear (dtest);
		for (i= 0; i < 4; i++)
		{
			for (j=0; j < 5; j++)
			{
				mpz_clear(res[i][j]);
				mpz_clear(res2_c[i][j]);
				mpz_clear(res2_d[i][j]);
			}
		}

	}

	mpz_mul(r_11,res[0][1],res[0][1]);
	mpz_mul(r_22,res[0][2],res[0][2]);
	mpz_mul(r_33,res[0][3],res[0][3]);
	mpz_mul(r_14,res[0][1],res[0][4]);
	mpz_mul(r_313,res[0][3],res[1][3]);
	mpz_mul(r_212,res[0][2],res[1][2]);
	mpz_mul(temp1,res[0][3],res[1][2]); 
	mpz_mul(temp2,res[0][2],res[1][3]);
	mpz_sub(diff,temp1,temp2);

/* Compute d0 */

	mpz_add(temp1,r_22,r_33); mpz_add(temp1,temp1,r_11); 
	mpz_mul(temp1,temp1,res[0][1]); 
	coef = -2;
	mpz_mul_si(d0,temp1,coef);

	if(!(*testr)) {

/* Compute d1 */

		mpz_add(temp1,r_313,r_212);
		coef = 2;
		mpz_mul_si(temp1,temp1,coef);
		mpz_sub(temp1,temp1,r_14);
		mpz_mul(d1,res[0][1],temp1);


/* Compute d2 */

		mpz_add(temp1,r_11,r_33);
		mpz_mul(temp1,temp1,res[1][2]);
		coef = -2;
		mpz_mul_si(temp1,temp1,coef);
		mpz_mul_si(temp2,r_313,coef);
		mpz_add(temp2,temp2,r_14);
		mpz_mul(temp2,temp2,res[0][2]);
		mpz_sub(d2,temp1,temp2);

/* Compute d3 */

		mpz_add(temp1,r_11,r_22);
		mpz_mul(temp1,temp1,res[1][3]);
		mpz_mul_si(temp1,temp1,coef);
		mpz_mul_si(temp2,r_212,coef);
		mpz_add(temp2,temp2,r_14);
		mpz_mul(temp2,temp2,res[0][3]);
		mpz_sub(d3,temp1,temp2);

/* Compute d4 */

		mpz_mul(temp1,res[0][3],res[3][4]);
		mpz_mul(temp2,res[0][2],res[2][4]);
		mpz_mul(temp3,res[0][1],res[1][4]);
		mpz_add(temp1,temp1,temp2); mpz_add(temp1,temp3,temp1);
		coef = 2;
		mpz_mul_si(temp1,temp1,coef);
		mpz_mul(temp1,temp1,res[0][1]);
		mpz_mul(temp2,res[1][3],res[1][3]);
		mpz_mul(temp3,res[1][2],res[1][2]);
		mpz_add(temp2,temp3,temp2);
		mpz_mul(temp2,temp2,res[0][1]);
		mpz_mul(temp3,res[2][3],diff);
		mpz_sub(temp2,temp3,temp2);
		coef = 4;
		mpz_mul_si(temp2,temp2,coef);
		mpz_add(d4,temp1,temp2);

/* Compute numerator of the radius of the smallest circumsphere of the edge */

		mpz_mul(temp1,d0,d4);
		mpz_mul(temp2,d3,d3);
		mpz_sub(temp2,temp2,temp1);
		mpz_mul(temp1,d2,d2);
		mpz_add(temp2,temp2,temp1);
		mpz_mul(temp1,d1,d1);
		mpz_add(num,temp1,temp2);

/* Compute denominator of the radius of the smallest circumsphere of the edge */

		mpz_mul(den,d0,d0);

/* check if radius is lower than ALPHA         */

		mpz_mul(temp1,den,alp);
		mpz_sub(temp2,num,temp1);

		if(mpz_sgn(temp2) < 0) (*testr)=1;
	}

/* Now check if edge (ab) is attached to c */

	if(!(*testa)) 
	{
		mpz_mul(temp1,res[1][2],res2_c[1][2]);
		mpz_mul(temp2,res[1][3],res2_c[1][3]);
		mpz_add(temp1,temp1,temp2);
		coef = -2;
		mpz_mul_si(temp1,temp1,coef);

		mpz_set_si(temp2,0);
		for (i=1; i<4; i++)
		{
			mpz_mul(temp3,res[0][i],res2_c[i][4]);
			mpz_add(temp2,temp2,temp3);
		}
		mpz_add(temp1,temp2,temp1);
		mpz_mul(temp1,temp1,res[0][1]);

		mpz_mul(temp2,res2_c[2][3],diff);
		mpz_mul_si(temp2,temp2,coef);
		mpz_sub(temp3,temp1,temp2);
		mpz_mul(dtest,temp3,d0);

		if(mpz_sgn(dtest) < 0) (*testa = 1);
	}

/* Now check if edge (ab) is attached to d */

	if(!(*testa)) 
	{ 
		mpz_mul(temp1,res[1][2],res2_d[1][2]);
		mpz_mul(temp2,res[1][3],res2_d[1][3]);
		mpz_add(temp1,temp1,temp2);
		coef = -2;
		mpz_mul_si(temp1,temp1,coef);

		mpz_set_si(temp2,0);
		for (i=1; i<4; i++)
		{
			mpz_mul(temp3,res[0][i],res2_d[i][4]);
			mpz_add(temp2,temp2,temp3);
		}
		mpz_add(temp1,temp2,temp1);
		mpz_mul(temp1,temp1,res[0][1]);

		mpz_mul(temp2,res2_d[2][3],diff);
		mpz_mul_si(temp2,temp2,coef);
		mpz_sub(temp3,temp1,temp2);
		mpz_mul(dtest,temp3,d0);

		if(mpz_sgn(dtest) < 0) (*testa) = 1;
	} 

/* 	Clear local GMP variables */

	mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3);
	mpz_clear (r_11); mpz_clear (r_22); mpz_clear (r_33); mpz_clear (r_14);
	mpz_clear (r_313); mpz_clear (r_212); mpz_clear (diff);
	mpz_clear (d0); mpz_clear (d1); mpz_clear (d2); mpz_clear (d3);
	mpz_clear (d4); 
	mpz_clear (num); mpz_clear (den); mpz_clear (dtest);

	for (i= 0; i < 4; i++)
	{
		for (j=0; j < 5; j++)
		{
			mpz_clear(res[i][j]);
			mpz_clear(res2_c[i][j]);
			mpz_clear(res2_d[i][j]);
		}
	}

}

/*------------------------------------------------------------------------*/
/* check_triangle_gmp:
   This program checks if a face of a tetrahedron belongs to the Alpha
   complex (with Alpha set to alp) and if this face is attached to the fourth
   vertex of the tetrahedron

	Input:

	For the three points m,n,p that form the triangles, the program
	needs as input the following determinants:

	S(i+j-2) = Minor(m,n,p,i,j,0)= det | m(i)  m(j)  1 |
					   | n(i)  n(j)  1 |
					   | p(i)  p(j)  1 |

	T(i) = Minor(m,n,p,i,4,0) = det | m(i)  m(4)  1 |
					| n(i)  n(4)  1 |
					| p(i)  p(4)  1 |

	U(i) = Minor(m,n,p,i,j,4) = det | m(i) m(j) m(4) |
					| n(i) n(j) n(4) |
					| p(i) p(j) p(4) |

If q is the fourth point of the tetrahedron,

	De1 = Minor(m,n,p,q,2,3,4,0)
	De2 = Minor(m,n,p,q,1,3,4,0)
	De3 = Minor(m,n,p,q,1,2,4,0)
	Dmnpq= Minor(m,n,p,q,1,2,3,0)
	Dmnp  = Minor(m,n,p,1,2,3)

	Output:

	testr	: flag set to 1 if ALPHA is larger than rho, the radius
		  of the circumsphere of the triangle
	testa	: flag set to 1 if the fourth point d of the tetrahedron
		  that contains the triangle (a,b,c) in inside the
		  circumsphere of (a,b,c)

*/

void check_triangle_gmp(mpz_t *S, mpz_t *T, mpz_t *U,
		mpz_t De1, mpz_t De2, mpz_t De3, mpz_t Dmnpq, mpz_t Dmnp,
		int *testr, int *testa, mpz_t alp)

{
	int i,coef;


/*	Local GMP variables */

	mpz_t temp1, temp2, temp3;
	mpz_t d0,d1,d2,d3,d4;
	mpz_t num,den,dtest;

/* 	Initialise local GMP variables */

	mpz_init(temp1); mpz_init(temp2); mpz_init(temp3);
	mpz_init (d0); mpz_init (d1); mpz_init (d2); mpz_init (d3);
	mpz_init (d4); 
	mpz_init (num); mpz_init (den); mpz_init (dtest);

/* First check if triangle is attached                               */

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

		if(mpz_sgn(dtest) > 0) (*testa)=1;
	}

	if(!(*testr)) {

/* Compute d0 */

		coef = 4;
		mpz_mul_si(d0,temp1,coef);

/* Compute d1 */

		mpz_mul(temp1,Dmnp,S[3]);
		coef = -2;
		mpz_mul_si(temp1,temp1,coef);
		mpz_mul(temp2,S[1],T[2]);
		mpz_add(temp1,temp2,temp1);
		mpz_mul(temp2,S[2],T[3]);
		mpz_add(temp1,temp2,temp1);
		mpz_mul_si(d1,temp1,coef);

/* Compute d2 */

		coef = 2;
		mpz_mul(temp1,Dmnp,S[2]);
		mpz_mul_si(temp1,temp1,coef);
		mpz_mul(temp2,S[3],T[3]);
		mpz_add(temp1,temp2,temp1);
		mpz_mul(temp2,S[1],T[1]);
		mpz_sub(temp1,temp2,temp1);
		mpz_mul_si(d2,temp1,coef);

/* Compute d3 */

		mpz_mul(temp1,Dmnp,S[1]);
		mpz_mul_si(temp1,temp1,coef);
		mpz_mul(temp2,S[2],T[1]);
		mpz_add(temp1,temp2,temp1);
		mpz_mul(temp2,S[3],T[2]);
		mpz_add(temp1,temp2,temp1);
		mpz_mul_si(d3,temp1,coef);

/* Compute d4 */

		mpz_mul(temp1,Dmnp,Dmnp);
		coef = -2;
		mpz_mul_si(temp1,temp1,coef);

		for (i=1; i<4; i++)
		{
			mpz_mul(temp2,S[i],U[i]);
			mpz_add(temp1,temp1,temp2);
		}
		coef = -4;
		mpz_mul_si(d4,temp1,coef);

/* Now compute numerator of the radius of the circumsphere of the triangle */

		mpz_mul(temp1,d0,d4);
		mpz_mul(temp2,d3,d3);
		mpz_sub(temp2,temp2,temp1);
		mpz_mul(temp1,d2,d2);
		mpz_add(temp2,temp2,temp1);
		mpz_mul(temp1,d1,d1);
		mpz_add(num,temp1,temp2);

/* Now compute denominator of the radius of the circumsphere of the triangle */

		mpz_mul(den,d0,d0);

/* Check if radius is lower than ALPHA */

		mpz_mul(temp1,den,alp);
		mpz_sub(temp2,num,temp1);

		if(mpz_sgn(temp2) < 0) (*testr) = 1;
	}

/* 	Clear local GMP variables */

	mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3);
	mpz_clear (d0); mpz_clear (d1); mpz_clear (d2); mpz_clear (d3);
	mpz_clear (d4); 
	mpz_clear (num); mpz_clear (den); mpz_clear (dtest);
}
/*------------------------------------------------------------------------*/
/* alf_tetra_gmp_	Version 1 11/24/2000	Patrice Koehl
 
 	This subroutine computes the radius R of the circumsphere containing
 	a tetrahedron [A,B,C,D], as well as check if any fourth point L 
 	(in A, B, C, D) of the tetrahedron is "hidden" by its opposite 
 	face [I,J,K], i.e. is interior to the cicumsphere of [I,J,K]
 
 	Since we are only interested at how R compares to Alpha, we don't
 	output R, rather the result of the comparison
 
 	Computation extends to all four faces of the tetrehedron, as
 	well as to all six edges of the tetrahedron
 
	This procedure works with Multiple Precision Integer Arithmetics  (MPIA)
 	The package GMP is used for MPIA (with a C wrapper)
*/ 
void alf_tetra_gmp_(double *a, double *b, double *c, double *d,
		double *ra, double *rb, double *rc, double *rd,
		int *testr, int *testa, double *scale, double *alpha)

{
	int i,j,k,coef,test_r,test_a,second[4];
	int ivalue;
	double value;

/*	Local GMP variables	*/

	mpz_t wa,wb,wc,wd;
	mpz_t num, den;
	mpz_t temp1,temp2,temp3;
	mpz_t val1,val2,val3,val4;
	mpz_t Dabc,Dabd,Dacd,Dbcd;
	mpz_t Det1,Det2,Det3,Det4,Dabcd;
	mpz_t alp;

	mpz_t a_mp[4], b_mp[4],c_mp[4],d_mp[4];
	mpz_t Sab[4],Sac[4],Sad[4],Sbc[4],Sbd[4],Scd[4];
	mpz_t Dab[4],Dac[4],Dad[4],Dbc[4],Dbd[4],Dcd[4];
	mpz_t Sa[4],Sb[4],Sc[4],Sd[4];
	mpz_t Sam1[4],Sbm1[4],Scm1[4],Sdm1[4];
	mpz_t Ta[4],Tb[4],Tc[4],Td[4];
	mpz_t Tam1[4],Tbm1[4],Tcm1[4],Tdm1[4];
	mpz_t Ua[4],Ub[4],Uc[4],Ud[4];
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

	mpz_init (num); mpz_init (den);

	mpz_init(temp1); mpz_init(temp2); mpz_init(temp3);
	mpz_init(val1); mpz_init(val2); mpz_init(val3); mpz_init(val4);

	for (i=0; i < 4; i++)
	{
		mpz_init(Sab[i]); mpz_init(Sac[i]); mpz_init(Sad[i]);
		mpz_init(Sbc[i]); mpz_init(Sbd[i]); mpz_init(Scd[i]);
		mpz_init(Dab[i]); mpz_init(Dac[i]); mpz_init(Dad[i]);
		mpz_init(Dbc[i]); mpz_init(Dbd[i]); mpz_init(Dcd[i]);
		mpz_init(Sa[i]); mpz_init(Sb[i]); mpz_init(Sc[i]); 
		mpz_init(Sd[i]);
		mpz_init(Sam1[i]); mpz_init(Sbm1[i]); mpz_init(Scm1[i]); 
		mpz_init(Sdm1[i]);
		mpz_init(Ta[i]); mpz_init(Tb[i]); mpz_init(Tc[i]); 
		mpz_init(Td[i]);
		mpz_init(Tam1[i]); mpz_init(Tbm1[i]); mpz_init(Tcm1[i]); 
		mpz_init(Tdm1[i]);
		mpz_init(Ua[i]); mpz_init(Ub[i]); mpz_init(Uc[i]); 
		mpz_init(Ud[i]);
		mpz_init(Deter[i]);
	}

	mpz_init(Dabc); mpz_init(Dabd); mpz_init(Dacd); mpz_init(Dbcd);
	mpz_init(Det1); mpz_init(Det2); mpz_init(Det3); mpz_init(Dabcd);
	mpz_init(Det4);

	mpz_init(alp);

/*	Transfer data in multiple precision */

	transfer_coord_to_mp(a,b,c,d,ra,rb,rc,rd,scale,a_mp,wa,
	b_mp,wb,c_mp,wc,d_mp,wd);
	value = (*alpha)*(*scale); ivalue = (int) floor(value); 
	mpz_set_si(alp,ivalue);

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
		mpz_neg(Sam1[i],Sa[i]); mpz_neg(Sbm1[i],Sb[i]);
		mpz_neg(Scm1[i],Sc[i]); mpz_neg(Sdm1[i],Sd[i]);
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
	We also need :
		Det = Det | m(1) m(2) m(3) m(4) |
			  | n(1) n(2) n(3) n(4) |
			  | p(1) p(2) p(3) p(4) |
			  | q(1) q(2) q(3) q(4) |
								*/

	mpz_mul(temp1,wa,Dbcd); mpz_mul(temp2,wb,Dacd);
	mpz_sub(temp3,temp2,temp1);
	mpz_mul(temp1,wc,Dabd); mpz_mul(temp2,wd,Dabc);
	mpz_sub(temp1,temp2,temp1); mpz_add(Dabcd,temp3,temp1);

/*
	The radius of the circumsphere of the weighted tetrahedron is then:
	r_t = (Det1*Det1 + Det2*Det2 + Det3*Det3 + 4*Det4*Dabcd)/(4*Det4*Det4)
								*/

	mpz_mul(temp1,Det4,Det4); coef=4; mpz_mul_si(den,temp1,coef);

	mpz_mul(temp1,Det1,Det1); mpz_mul(temp2,Det2,Det2);
	mpz_add(temp1,temp1,temp2); mpz_mul(temp2,Det3,Det3);
	mpz_add(temp1,temp1,temp2); mpz_mul(temp2,Det4,Dabcd);
	mpz_mul_si(temp2,temp2,coef); mpz_add(num,temp1,temp2);
	
	mpz_mul(temp1,den,alp); mpz_sub(temp2,num,temp1);
/* 
 	If tetrahedron is part of the alpha shape, then the 4 triangles,
 	the 6 edges and the four vertices are also part of the alpha
 	complex
 								*/
	if(mpz_sgn(temp2) < 0) 
	{
		for (i=0; i<15;i++) 
		{
			testr[i]=1;
			testa[i]=1;
		}
	}
	else
	{
/* 
 	Now check all four faces of the tetrahedra:
 	We look both if the fourth vertex is "hidden" by the face of 
 	interest, and then we compute the radius of the circumsphere
 	of the face
 
 	We first need the minors:
 
 	Tm(i) = Minor(n,p,q,i,4,0) = det | n(i)  n(4)  1 |
 					 | p(i)  p(4)  1 |
 					 | q(i)  q(4)  1 |
 
 	Um(i) = Minor(n,p,q,i,j,4) = det | n(i) n(j) n(4) |
 					 | p(i) p(j) p(4) |
 					 | q(i) q(j) q(4) |
 									*/

		for (i=1; i<4; i++)
		{
			mpz_sub(Dab[i],a_mp[i],b_mp[i]);
			mpz_sub(Dac[i],a_mp[i],c_mp[i]);
			mpz_sub(Dad[i],a_mp[i],d_mp[i]);
			mpz_sub(Dbc[i],b_mp[i],c_mp[i]);
			mpz_sub(Dbd[i],b_mp[i],d_mp[i]);
			mpz_sub(Dcd[i],c_mp[i],d_mp[i]);
		}
		for (i=1; i<4; i++)
		{
			mpz_mul(temp1,wb,Dcd[i]); mpz_mul(temp2,wd,Dbc[i]);
			mpz_mul(temp3,wc,Dbd[i]); mpz_add(temp1,temp1,temp2);
			mpz_sub(Ta[i],temp3,temp1); 
			mpz_neg(Tam1[i],Ta[i]);
			mpz_mul(temp1,wa,Dcd[i]); mpz_mul(temp2,wd,Dac[i]);
			mpz_mul(temp3,wc,Dad[i]); mpz_add(temp1,temp1,temp2);
			mpz_sub(Tb[i],temp3,temp1); 
			mpz_neg(Tbm1[i],Tb[i]);
			mpz_mul(temp1,wa,Dbd[i]); mpz_mul(temp2,wd,Dab[i]);
			mpz_mul(temp3,wb,Dad[i]); mpz_add(temp1,temp1,temp2);
			mpz_sub(Tc[i],temp3,temp1); 
			mpz_neg(Tcm1[i],Tc[i]);
			mpz_mul(temp1,wa,Dbc[i]); mpz_mul(temp2,wc,Dab[i]);
			mpz_mul(temp3,wb,Dac[i]); mpz_add(temp1,temp1,temp2);
			mpz_sub(Td[i],temp3,temp1); 
			mpz_neg(Tdm1[i],Td[i]);
			mpz_mul(temp1,wb,Scd[i]); mpz_mul(temp2,wd,Sbc[i]);
			mpz_mul(temp3,wc,Sbd[i]); mpz_add(temp1,temp1,temp2);
			mpz_sub(Ua[i],temp1,temp3);
			mpz_mul(temp1,wa,Scd[i]); mpz_mul(temp2,wd,Sac[i]);
			mpz_mul(temp3,wc,Sad[i]); mpz_add(temp1,temp1,temp2);
			mpz_sub(Ub[i],temp1,temp3);
			mpz_mul(temp1,wa,Sbd[i]); mpz_mul(temp2,wd,Sab[i]);
			mpz_mul(temp3,wb,Sad[i]); mpz_add(temp1,temp1,temp2);
			mpz_sub(Uc[i],temp1,temp3);
			mpz_mul(temp1,wa,Sbc[i]); mpz_mul(temp2,wc,Sab[i]);
			mpz_mul(temp3,wb,Sac[i]); mpz_add(temp1,temp1,temp2);
			mpz_sub(Ud[i],temp1,temp3);
		}
		mpz_neg(val1,Det1); mpz_neg(val2,Det2);
		mpz_neg(val3,Det3); mpz_neg(val4,Det4);

/* First check face abc */

		if(!testa[1])
		{
			test_r = testr[1]; test_a = testa[1];
			second[0]=0; if(test_r) second[0]=1;
			check_triangle_gmp(Sd,Td,Ud,Det1,Det2,Det3,Det4,Dabc,
				&test_r,&test_a,alp);
			testr[1]=test_r; if(!testa[1]) testa[1]=test_a;
			if(testa[1]) testr[1]=0;
		}

/* Now check face abd */
 
		if(!testa[2])
		{
			test_r = testr[2]; test_a = testa[2];
			second[1] = 0; if(test_r) second[1]=1;
			check_triangle_gmp(Sc,Tc,Uc,val1,val2,val3,val4,Dabd,
				&test_r,&test_a,alp);
			testr[2]=test_r; if(!testa[2]) testa[2]=test_a;
			if(testa[2]) testr[2]=0;
		}

/* Now check face acd */
 
		if(!testa[3])
		{
			test_r = testr[3]; test_a = testa[3];
			second[2] = 0; if(test_r) second[2]=1;
			check_triangle_gmp(Sb,Tb,Ub,Det1,Det2,Det3,Det4,Dacd,
				&test_r,&test_a,alp);
			testr[3]=test_r; if(!testa[3]) testa[3]=test_a;
			if(testa[3]) testr[3]=0;
		}
 
/* Now check face bcd */
 
		if(!testa[4])
		{
			test_r = testr[4]; test_a = testa[4];
			second[3] = 0; if(test_r) second[3]=1;
			check_triangle_gmp(Sa,Ta,Ua,val1,val2,val3,val4,Dbcd,
				&test_r,&test_a,alp);
			testr[4]=test_r; if(!testa[4]) testa[4]=test_a;
			if(testa[4]) testr[4]=0;
		}
/* 
 	Now consider edges:
 	Start by checking each triangle: if it belongs to the alpha shape,
 	so does its 3 edges:
									*/
		if((testr[1]!=0) && (second[0]!=0) ) {
			testr[5]=1;testr[6]=1;testr[8]=1;
			testa[5]=1;testa[6]=1;testa[8]=1;
		}
		if((testr[2]!=0) && (second[1]!=0) ) {
			testr[5]=1;testr[7]=1;testr[9]=1;
			testa[5]=1;testa[7]=1;testa[9]=1;
		}
		if((testr[3]!=0) && (second[2]!=0) ) {
			testr[6]=1;testr[7]=1;testr[10]=1;
			testa[6]=1;testa[7]=1;testa[10]=1;
		}
		if((testr[4]!=0) && (second[3]!=0) ) {
			testr[8]=1;testr[9]=1;testr[10]=1;
			testa[8]=1;testa[9]=1;testa[10]=1;
		}
 
/* Now check each edges */
 
		if(!testa[5])
		{
			test_r = testr[5]; test_a = testa[5];
			check_edge_gmp(a_mp,b_mp,wa,wb,Dab,Sab,Sd,Sc,
			Td,Tc,&test_r,&test_a,alp); 
			testr[5]=test_r; if(!testa[5]) testa[5]=test_a;
			if(testa[5]) testr[5]=0;
		}
		if(!testa[6])
		{
			test_r = testr[6]; test_a = testa[6];
			check_edge_gmp(a_mp,c_mp,wa,wc,Dac,Sac,Sdm1,Sb,
			Tdm1,Tb,&test_r,&test_a,alp); 
			testr[6]=test_r; if(!testa[6]) testa[6]=test_a;
			if(testa[6]) testr[6]=0;
		}
		if(!testa[7])
		{
			test_r = testr[7]; test_a = testa[7];
			check_edge_gmp(a_mp,d_mp,wa,wd,Dad,Sad,Scm1,Sbm1,
			Tcm1,Tbm1,&test_r,&test_a,alp); 
			testr[7]=test_r; if(!testa[7]) testa[7]=test_a;
			if(testa[7]) testr[7]=0;
		}
		if(!testa[8])
		{
			test_r = testr[8]; test_a = testa[8];
			check_edge_gmp(b_mp,c_mp,wb,wc,Dbc,Sbc,Sd,Sa,
			Td,Ta,&test_r,&test_a,alp); 
			testr[8]=test_r; if(!testa[8]) testa[8]=test_a;
			if(testa[8]) testr[8]=0;
		}
		if(!testa[9])
		{
			test_r = testr[9]; test_a = testa[9];
			check_edge_gmp(b_mp, d_mp, wb, wd,Dbd,Sbd,Sc,Sam1,
			Tc,Tam1,&test_r,&test_a,alp); 
			testr[9]=test_r; if(!testa[9]) testa[9]=test_a;
			if(testa[9]) testr[9]=0;
		}
		if(!testa[10])
		{
			test_r = testr[10]; test_a = testa[10];
			check_edge_gmp(c_mp,d_mp,wc,wd,Dcd,Scd,Sb,Sa,
			Tb,Ta,&test_r,&test_a,alp); 
			testr[10]=test_r; if(!testa[10]) testa[10]=test_a;
			if(testa[10]) testr[10]=0;
		}

/*	Set all vertices in alpha complex */

		for (i=11;i<15;i++) testr[i]=1;
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

	mpz_clear (num); mpz_clear (den);

	mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3);

	for (i=0; i < 4; i++)
	{
		mpz_clear(Sab[i]); mpz_clear(Sac[i]); mpz_clear(Sad[i]);
		mpz_clear(Sbc[i]); mpz_clear(Sbd[i]); mpz_clear(Scd[i]);
		mpz_clear(Dab[i]); mpz_clear(Dac[i]); mpz_clear(Dad[i]);
		mpz_clear(Dbc[i]); mpz_clear(Dbd[i]); mpz_clear(Dcd[i]);
		mpz_clear(Sa[i]); mpz_clear(Sb[i]); mpz_clear(Sc[i]); 
		mpz_clear(Sd[i]);
		mpz_clear(Sam1[i]); mpz_clear(Sbm1[i]); mpz_clear(Scm1[i]); 
		mpz_clear(Sdm1[i]);
		mpz_clear(Ta[i]); mpz_clear(Tb[i]); mpz_clear(Tc[i]); 
		mpz_clear(Td[i]);
		mpz_clear(Tam1[i]); mpz_clear(Tbm1[i]); mpz_clear(Tcm1[i]); 
		mpz_clear(Tdm1[i]);
		mpz_clear(Ua[i]); mpz_clear(Ub[i]); mpz_clear(Uc[i]); 
		mpz_clear(Ud[i]);
		mpz_clear(Deter[i]);
	}

	mpz_clear(Dabc); mpz_clear(Dabd); mpz_clear(Dacd); mpz_clear(Dbcd);
	mpz_clear(Det1); mpz_clear(Det2); mpz_clear(Det3); mpz_clear(Dabcd);
	mpz_clear(Det4);

	mpz_clear(alp);

}
