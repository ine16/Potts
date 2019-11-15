#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pcg_basic.h"

/*Modelo de Potts usando heat bath (codigo no muy compacto pero bastante legible)*/


void inicio(int *centinela, int size, int num_estados){
	//pcg32_random_t rngptr;
	for(int cont=0;cont<size;cont++){
		//*centinela=pcg32_boundedrand_r(&rngptr, num_estados)+1;
		*centinela=1;
		centinela++;
	}
}

void bolsman(double b, double *vect1){
	*vect1=1;
	for (int ind=1;ind<5;ind++){
		*(vect1+ind)=exp(b*ind);
	}
}

int kronecker(int sitio1, int sitio2){
	if (sitio1==sitio2)
		return 1;
	else
		return 0;
}

void prob_acumulativa(double *expo, double *vect2, int estados, int izq, int der, int arriba, int abajo){
	//Primero calcula los factores de boltzmann
	for (int cc=1; cc<(estados+1); cc++){
		*(vect2+cc-1)= *(expo+(kronecker(cc,izq)+kronecker(cc,der)+kronecker(cc,arriba)+kronecker(cc,abajo)));
	}
	double suma_boltz=0;
	//Ahora la traza total (fx particion)
	for (int cc1=0; cc1<estados; cc1++){
		suma_boltz += *(vect2+cc1);
	}
	//Despues arma las probabilidades de cada estado dividiendo por la traza
	for (int cc2=0; cc2<estados; cc2++){
		*(vect2+cc2)= *(vect2+cc2)/suma_boltz;
	}
	//Y por ultimo las suma para armar la probabilidad acumulativa
	for (int cc3=1; cc3<estados; cc3++){
		*(vect2+cc3)= *(vect2+cc3)+*(vect2+cc3-1);
	}
	*(vect2+estados)=1;
}



int mod (int a, int b) //robe esta implementacion de mod porque C hace cosas locas con los negativos
{
   if (b < 0) //you can check for b == 0 separately and do what you want
     return mod(a, -b);   
   int ret = a % b;
   if (ret < 0)
     ret+=b;
   return ret;
}

double maspoblado(double *proporciones, int tamanio){
	double bicho = 0;
	for(int cont1=0; cont1<tamanio; cont1++){
		if ((*(proporciones+cont1))>bicho)
			bicho = *(proporciones+cont1);
	}
	return bicho;
}

double bootstrap(int *datos, int total_datos, double inv_T, int sitios_red){
	pcg32_random_t rngpunt;
	int elemento, items, n;
	double mean, mean_sq, aux, mean_items, mean_sq_items, sigma;
	mean=0;
	mean_items=0;
	mean_sq=0;
	mean_sq_items=0;
	n=50;//numero de veces que agarro subconjuntos para calcular
	items=300;//subconjuntos que uso para calcular los
	for (int contador1=0; contador1<n;contador1++){
		for (int contador2=0; contador2<items; contador2++){
			elemento = pcg32_boundedrand_r(&rngpunt, total_datos);
			mean_items += *(datos+elemento);
			mean_sq_items += pow(*(datos+elemento),2.0);
		}
		mean_items = mean_items / items;
		mean_sq_items = mean_sq_items / items;
		aux = inv_T*(mean_sq_items - pow(mean_items, 2.0))/sitios_red;
		mean += aux;
		mean_sq += pow(aux,2.0);
	}
	mean = mean/n;
	mean_sq = mean_sq/n;
	sigma = sqrt(mean_sq - pow(mean,2.0));
	return sigma;
}

void iniciar(int *vect3, int length){
	for(int val=0; val<length; val++){
		*vect3 = 0;
		vect3++;
	}
}

double Magnetizacion(int *ptr_magn, double *salida_magn, int sitiosm){
	double m=0;
	for (int pituto=0; pituto<sitiosm; pituto++){
		m += *(ptr_magn+pituto);
	}
	*salida_magn += m;
	*(salida_magn+1) += pow(m,2);
	return pow(m,2);
}


int main(){
	//Hamiltoniano H=-J*sum(delta_kron(i,j)), tomo J=1
	//Como E<=0 y las proba son prop a exp(-beta*E) voy a volar los signos (-) de H y de la exp
	//Cuando calcule la energia ahi si voy a agregar el signo (-)

	int q=2; //numero de estados posibles para cada sitio (valor minimo: 2) TOMAR 2,4,10
	int L=20; //"largo" del lado de la red cuadrada PARA HACER POSTA PONER L=50,30,80
	int N=L*L; //numero de sitios
	int pasos_termalizacion=5*N; //numero de pasos que pienso dejar pasar antes de tomar promedios (antes decia 2N)
	int MCS=10*N; //numero de pasos Montecarlo para el programa (antes decia 6*N)
	int medidas = MCS - pasos_termalizacion; //por lo que divido para calcular el promedio

	/* declaro la red */
	int red[L][L];
	/* inicializo la red aleatoriamente, asignando a cada sitio un valor entre 1,...,q */
	inicio(&red[0][0],N,q);

	pcg32_random_t rngptr1, rngptr2, myrng; //cosos para el generador de aleatorios
	/*COMO LLAMAR AL GENERADOR DE ALEATORIOS PARA QUE FUNCIONE:
	//pcg32_random_t rng1, rngptr, myrng; //esto va si o si pa que arranque el generador pcg
	//pcg32_srandom_r(&rng1, time(NULL), (intptr_t)&rng1); //tira num random entero con andasaber que limites
	//int i = pcg32_boundedrand_r(&rngptr, cut); //tira num random entero q, 0<=q<cut
	//double valor_rand = ldexp(pcg32_random_r(&myrng), -32); //tira random tipo double entre 0 y 1
	*/

	int i,j; //indices que van a indicar el sitio (i,j) a modificar
	int newstate; //variable que voy a usar para definir nuevo estado del sitio (i,j)
	FILE *pointer_medidas; //el puntero del archivo que tiene todas las medidas
	FILE *pointer_medidasE; //el puntero del archivo que tiene todas las medidas




	/* --- ACA ARRANCA EL MONTECARLO COMPLETO --- */
	pointer_medidas= fopen("medidas_q2_L20.txt", "w");
	pointer_medidasE= fopen("datosE_q2_L20.txt", "w");
	int datos=301;
	for(int xx=1; xx<datos; xx++){
		printf("%d\n",xx);
		double beta = 5 - xx * 4.75 / datos ; //beta=1./(k*T) con k=1 //el original era con 4.95

		/* Variables a medir para cada valor de beta */
		double E=0; //energia
		double cv=0; //calor especifico
		double M[2]={0,0}; //M(0) es la magnetizacion, M(1) es la susceptibilidad
		int M_vec[medidas], subtotales[q];
		double chi_err;
		double EEE;
		
		int max_fraction=0;
		double order_param;

		double exp[5];
		bolsman(beta, exp); //dejo calculados los valores de exp que voy a usar

		//inicializo en cero a los vectores para que guarden los valores para el calculo de errores
		iniciar(M_vec,medidas);

		/* --- ACA ARRANCARIA EL MONTECARLO PARA CADA VALOR DE BETA--- */
		for(int pasos=0; pasos<MCS; pasos++){
			for(int cuchuflo=0; cuchuflo<N; cuchuflo++){
				/* elijo sitio (i,j) */
				int i = pcg32_boundedrand_r(&rngptr1, L);
				int j = pcg32_boundedrand_r(&rngptr2, L);
				//printf("%d ",i);
				double valor_rand = ldexp(pcg32_random_r(&myrng), -32);
				/* HEAT BATH */
				double proba[q+1];
				prob_acumulativa(exp,proba,q,red[i][mod(j-1,L)],red[i][mod(j+1,L)],red[mod(i-1,L)][j],red[mod(i+1,L)][j]);
				newstate=0;
				while(proba[newstate]<valor_rand){
					newstate++;}
				red[i][j]=newstate+1;
				///---------------
				/* METROPOLIS
				lugar=0;
				for(int pepo=1; pepo<(q+1); pepo++){
					if (red[i][j] != pepo){
						estados_posibles[lugar]=pepo;
						lugar++;
					}
				}
				newstate= estados_posibles[pcg32_boundedrand_r(&rngptr3, q)+1];
				deltaE=-(kronecker(newstate,red[i][mod(j-1,L)])+kronecker(newstate,red[i][mod(j+1,L)])+kronecker(newstate,red[mod(i-1,L)][j])+kronecker(newstate,red[mod(i+1,L)][j]))+(kronecker(red[i][j],red[i][mod(j-1,L)])+kronecker(red[i][j],red[i][mod(j+1,L)])+kronecker(red[i][j],red[mod(i-1,L)][j])+kronecker(red[i][j],red[mod(i+1,L)][j]));
				if (deltaE<0 || valor_rand<exp_M[deltaE+4])
					red[i][j]= newstate;*/
				///---------------
			}////ACA TERMINA UN PASO MONTECARLO

			/* Medidas de variables termodinamicas y parametro de orden */
			
			/*
			for (int ii=0;ii<L;ii++){
				for (int jj=0;jj<L;jj++){
					EEE = EEE - kronecker(red[ii][jj],red[ii][mod(jj-1,L)]) - kronecker(red[ii][jj],red[mod(ii-1,L)][jj]);
				}
			}*/
			if (pasos >= (pasos_termalizacion - 1)){
				iniciar(subtotales,q);
				EEE=0; //variable auxiliar para poder calcular E y cv
				for (int ii=0;ii<L;ii++){
					for (int jj=0;jj<L;jj++){
						EEE = EEE - kronecker(red[ii][jj],red[ii][mod(jj-1,L)]) - kronecker(red[ii][jj],red[mod(ii-1,L)][jj]);
						subtotales[red[ii][jj]-1] += 1;
					}
				}
				E += EEE;
				fprintf(pointer_medidasE, "%f ", EEE);
				max_fraction += maspoblado(&subtotales[0],q);
				cv += pow(EEE,2);
				M_vec[pasos-pasos_termalizacion+1] = Magnetizacion(&red[0][0], M, N);
			}
		}//ACA TERMINA EL MONTECARLO PARA CADA VALOR DE BETA
	
		//EN DATOS FINALES DE L=50, DIVIDIR TODAS LAS MEDIDAS POR N Y PONER T EN VEZ DE BETA
		E = (double) E/medidas;
		cv = pow(beta,2)*(cv/medidas - pow(E, 2));
		M[0] = M[0]/medidas;
		M[1] = beta*(M[1]/medidas - pow(M[0], 2));
		chi_err = bootstrap(M_vec,medidas,beta,N);
		order_param = (double) (q * max_fraction/medidas/N - 1)/(q-1);
		//----ESTA PARTE CAMBIA RESPECTO DEL DE L=50
		E = E/N;
		cv = cv/N;
		M[0] = M[0]/N;
		M[1] = M[1]/N;
		//----
		double T=1./beta;
		fprintf(pointer_medidas, "%f ", T);
		fprintf(pointer_medidas, "%f ", E);
		fprintf(pointer_medidas, "%f ", cv);
		fprintf(pointer_medidas, "%f ", M[0]);
		fprintf(pointer_medidas, "%f ", M[1]);
		fprintf(pointer_medidas, "%f ", chi_err);
		fprintf(pointer_medidas, "%f\n", order_param);
		fprintf(pointer_medidasE, "\n");
	}
	////ACA CERRARIA EL FOR DEL BETA
	fclose(pointer_medidas);
	fclose(pointer_medidasE);
	return 0;
}
