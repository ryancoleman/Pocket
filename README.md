		Detecting Pockets in Proteins

1. Background

Active sites of proteins are often found in structural pockets or 
cavities. The Alpha Shape theory (see http://www.alphashapes.org/) 
provides an analytical method for detecting pockets in proteins, and 
measuring their volume and surface. The method has been fully 
described in:
Jie Liang, Herbert Edelsbrunner, and Clare Woodward. 1998. Anatomy 
of Protein Pockets and Cavities: Measurement of Binding Site Geometry
and Implications for Ligand Design. Protein Science, 7:1884-1897

This method was implemented in a new program, called Pocket. 

2. The Program

We distribute the source code for the whole program, including all
subroutines, under the LGPL licensing code. Please read the text of
the license (provided with the distribution). Basically, you are
allowed to use and modify this code freely. You are also allowed to
distribute it, provided that you provide the source code for all
functions involded in Pocket

The Makefile provided with the distribution is intentionally kept
simple. The only requirements to compile and run Pocket is that you
have a version of the GNU Multiple Precision (GMP) library installed
on your system.

For simplicity, we also include in the distribution one executable,
pocket, compiled statically under Linux 8.0 using the Fortran and C Intel
compilers for Linux.

Before (re)compiling the program, please check the file pocket.h:
it provides upper limits for the number of points as well as the
number of simplices the program can handle. Please modify this
file based on your own needs.

3. Running the Program

Pocket is a stand-alone Fortran program. To run Pocket, just type:
./pocket
The program then asks you 4 questions:

	- name of the PDB file describing the protein structure
	- name of the file defining atomic radii
	- generic name for output files
	- radius of the probe

(a) The PDB file should be in standard PDB format. Hydrogens are
excluded from the calculation

(b) Standard atomic radii for all types of atoms are provided in
the distribution in a file called opls_rad.dat . Please edit this
file to your own taste

(c) Based on the generic name (for example "test") you provide to
the program, it will generate 6 different output files:

	test.info	contains a summary of all calculations:
			surface and volume of the protein, number
			of pockets and cavities found, their
			surface and volume, as well as information
			on their mouths

	test.atvol	For all atoms of the protein, this file
			provide their contributions to the 
			total surface and volume of the protein,
			We also provide the derivatives of the
			total surface and total volume with 
			repsect to the coordinates of each atom

	test.pocket	list all atoms bordering the pockets of
			the protein

	test.mouth	list all atoms defining the mouths of
			the pockets

	test.rasmol	a rasmol script allowing to visualise
			the pockets. If you have Rasmol installed,
			use this script by typing:
			rasmol -script test.rasmol test.pdb

	test.kin	a Kinemage script that allows you to
			visualise different eometric properties
			of you protein. You need the program
			"mage".

4. A test case

We provide with the distribution a test PDB file called small.pdb,
as well as its corresponding output files with generic name small.

5. Disclaimer

The program is provided "as is". Contact koehl@csb.stanford.edu
if you have any problems.

6. Changes by Ryan Coleman

Writes out alpha shape data structure for use/reading in by other code.
See small.alpha for sample file

[![doi](https://zenodo.org/badge/3853/ryancoleman/Pocket.png)](http://dx.doi.org/10.5281/zenodo.10212)
