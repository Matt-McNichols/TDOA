#include<stdio.h>
#include<math.h>
#include "utilityFunctions.h"
#include<gsl/gsl_linalg.h>

void linSys(gsl_matrix *A, gsl_vector *B, gsl_vector *X ){
	int s;
	gsl_matrix *A_orriginal = A;
	gsl_vector *residual = gsl_vector_alloc(3);
	gsl_permutation *perm = gsl_permutation_alloc(3);
	gsl_linalg_LU_decomp(A,perm,&s);
	gsl_linalg_LU_refine(A_orriginal,A,perm, B, X, residual);
	return;
}
/*
void TEST_linSys(){
	// test solving linear system of equations
	double a_data[] = { 0.18, 0.60, 0.57, 0.96,
                      0.41, 0.24, 0.99, 0.58,
                      0.14, 0.30, 0.97, 0.66,
                      0.51, 0.13, 0.19, 0.85 };

	double b_data[] = { 1.0, 2.0, 3.0, 4.0 };

  	gsl_matrix_view a = gsl_matrix_view_array (a_data, 4, 4);

	gsl_vector_view b = gsl_vector_view_array (b_data, 4);

	gsl_vector *x = gsl_vector_alloc (4);

	linSys(a,b,x);
	printf("X = \n");
	gsl_vector_fprintf(stdout, x,"%g");
}
*/
void calcR(gsl_matrix_view P,gsl_vector_view E, gsl_vector *R){
	double e0 = gsl_vector_get(&E.vector,0);
	double e1 = gsl_vector_get(&E.vector,1);
	double e2 = gsl_vector_get(&E.vector,2);
	double r0 = sqrt(
			pow(gsl_matrix_get(&P.matrix,0,0) + e0,2) +
			pow(gsl_matrix_get(&P.matrix,0,1) + e1,2) +
			pow(gsl_matrix_get(&P.matrix,0,2) + e2,2) 
			);
	double r1 = sqrt(
			pow(gsl_matrix_get(&P.matrix,1,0) + e0,2) +
			pow(gsl_matrix_get(&P.matrix,1,1) + e1,2) +
			pow(gsl_matrix_get(&P.matrix,1,2) + e2,2) 
			);
	double r2 = sqrt(
			pow(gsl_matrix_get(&P.matrix,2,0) + e0,2) +
			pow(gsl_matrix_get(&P.matrix,2,1) + e1,2) +
			pow(gsl_matrix_get(&P.matrix,2,2) + e2,2) 
			);
	double r3 = sqrt(
			pow(gsl_matrix_get(&P.matrix,3,0) + e0,2) +
			pow(gsl_matrix_get(&P.matrix,3,1) + e1,2) +
			pow(gsl_matrix_get(&P.matrix,3,2) + e2,2) 
			);
	double r4 = sqrt(
			pow(gsl_matrix_get(&P.matrix,4,0) + e0,2) +
			pow(gsl_matrix_get(&P.matrix,4,1) + e1,2) +
			pow(gsl_matrix_get(&P.matrix,4,2) + e2,2) 
			);
	gsl_vector_set(R,0,r0);
	gsl_vector_set(R,1,r1);
	gsl_vector_set(R,2,r2);
	gsl_vector_set(R,3,r3);
	gsl_vector_set(R,4,r4);
}

void calcTAU(gsl_matrix_view P,gsl_vector *T, double v, gsl_vector *TAU){
	double To = gsl_vector_get(T,0);
	printf("To = %g\n",To);
	double val;
	int m;
	for(m = 1; m<5; m++){
		val = v*(gsl_vector_get(T,m)-To);
		gsl_vector_set(TAU,m-1,val);
	}
}

void constCoef(gsl_matrix_view P,gsl_vector *TAU, gsl_matrix *COEF, gsl_vector *D){
	/*
	Am = (2*Xm/TAUm) - (2*X1/TAU1)
	bm = (2*Ym/TAUm) - (2*Y1/TAU1)
	cm = (2*Zm/TAUm) - (2*Z1/TAU1)
	dm = TAUm - TAU1 - (Xm^2 + Ym^2 +Zm^2)/TAUm + (X1^2 + Y1^2 +Z1^2)/TAU1
	*/
	double TAUo = gsl_vector_get(TAU,0);
	double P1x = gsl_matrix_get(&P.matrix,1,0);
	double P1y = gsl_matrix_get(&P.matrix,1,1);
	double P1z = gsl_matrix_get(&P.matrix,1,2);
	double subA = 2*P1x/TAUo;
	double subB = 2*P1y/TAUo;
	double subC = 2*P1z/TAUo;
	double subD = (pow(P1x,2) + pow(P1y,2) + pow(P1z,2))/TAUo;
	int m;
	double Am, Bm, Cm, Dm, Px, Py, Pz, tmp;
	for(m = 2; m<5; m++){
		Px = gsl_matrix_get(&P.matrix,m,0);
		Py = gsl_matrix_get(&P.matrix,m,1);		
		Pz = gsl_matrix_get(&P.matrix,m,2);
		Am = (2*Px/gsl_vector_get(TAU,m-1)) - subA;
		Bm = (2*Py/gsl_vector_get(TAU,m-1)) - subB;
		Cm = (2*Pz/gsl_vector_get(TAU,m-1)) - subC;
		tmp = (pow(Px,2) + pow(Py,2) + pow(Pz,2))/gsl_vector_get(TAU,m-1);
		Dm = gsl_vector_get(TAU,m-1) - TAUo - tmp +subD;
		gsl_matrix_set(COEF,m-2,0,Am);
		gsl_matrix_set(COEF,m-2,1,Bm);
		gsl_matrix_set(COEF,m-2,2,Cm);
		gsl_vector_set(D,m-2,Dm);m = (2*Xm/TAUm) - (2*X1/TAU1)
    # bm = (2*Ym/TAUm) - (2*Y1/TAU1)
    # cm = (2*Zm/TAUm) - (2*Z1/TAU1)
    # dm = TAUm - TAU1 - (Xm^2 + Ym^2 +Zm^2)/TAUm + (X1^2 + Y1^2 +Z1^2)/TAU1

		

	}


}


void problem1(gsl_matrix_view P, gsl_vector_view E, double v, gsl_vector *T){
	gsl_vector *R = gsl_vector_alloc(5);
	double vInv = 1/v;
	printf("calc R\n");
	calcR(P,E,R);
	gsl_vector_fprintf(stdout,R,"%g");
	*T = *R;
	gsl_vector_scale(T,vInv);
}

void problem2(gsl_matrix_view P, gsl_vector * T, double v, gsl_vector *approxE){
	gsl_vector *TAU = gsl_vector_alloc(4);
	gsl_vector *D   = gsl_vector_alloc(3);
	gsl_matrix *COEF = gsl_matrix_alloc(3,3);
	calcTAU(P,T,v,TAU);
	printf("TAU = \n");
	gsl_vector_fprintf(stdout,TAU,"%g");
	constCoef(P,TAU,COEF,D);
	printf("COEF = \n");
	gsl_matrix_fprintf(stdout,COEF,"%g");
	printf("D = \n");
	gsl_vector_fprintf(stdout,D,"%g");
	linSys(COEF, D, approxE);
}

int main(){
	printf("\n\nParameters\n--------------------\n");
	double v = 1481;
	double P_data [] = { 0, 0, 1,
			     1, 0, 0,
			     0, 1, 0,
			    -1, 0, 0,
			     0,-1, 0 };
	double E_data[] = { -4, 6, 5 };
	gsl_vector_view E = gsl_vector_view_array(E_data, 3);
	gsl_matrix_view P = gsl_matrix_view_array(P_data, 5, 3);
	gsl_vector *T = gsl_vector_alloc(5);
	gsl_vector *approxE = gsl_vector_alloc(3);
	//gsl_matrix_fprintf(stdout, &P.matrix, "%g");
	//printf("TEST linSys\n--------------------\n");
	//TEST_linSys();
	printf("\nProblem1: calc T\n--------------------\n");
	problem1(P,E,v,T);
	printf("T = \n");
	gsl_vector_fprintf(stdout,T,"%g");
	printf("\nProblem2: approximate location of E\n--------------------\n");
	problem2(P,T,v,approxE);
	//printf("approxE = \n");
	//gsl_vector_fprintf(stdout,approxE,"%g");
	
	 
	return 0;
}

