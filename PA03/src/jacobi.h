/* JACOBI.H

   NOTES: Header file to be included in every module, that needs
	  to call jacobi().

*/
typedef unsigned dimension;
int jacobi(double **S, dimension n, double *w, double **V);
void jacobi_set_max_iterations(unsigned long iter);

