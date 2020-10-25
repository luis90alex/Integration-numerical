#include "Myheader.h"

double a, b, u, e;

void trapezi(double fi0, double g, double L, double e, double u){
    double h_trapecio, periodo_armonico_t, error;
    double w0;
    int i , nt, j;

    FILE *archivo_trapecio;
	archivo_trapecio = fopen("calculo_trapecio.txt", "w");

	printf ("Este programa calcula la integral a traves del metodo de los trapecios.\n\n\n");

	printf ("Define el numero de divisiones n: ");///seria nuestro número máximo de iteraciones
	scanf ("%d", &nt);

	fprintf (archivo_trapecio, "PROGRAMA REALIZADO POR:\n   Adrià Medeiros (NIU:1388403)\n   Luís Farfan (NIU:1388330)\n   Lucas Fuentealba(NIU:1388326)\n\n");
	fprintf (archivo_trapecio, "CONDICIONES INICIALES:\n");
    fprintf (archivo_trapecio, "\n- Numero de divisiones n: %d", nt);
	fprintf (archivo_trapecio, "\n- Condicion inicial fi0: %.15lf", fi0);
	fprintf (archivo_trapecio, "\n- Gravedad g (m/s^2): %.15lf", g);
	fprintf (archivo_trapecio, "\n- Longitud L (m): %.15lf", L);
	fprintf (archivo_trapecio, "\n- Tolerancia e: %.15lf", e);

	w0 = sqrt(g/L);
    periodo_armonico_t=2*PI*sqrt(L/g);
	a = 0;
	b = PI/2;

	double *periodo;
	double *primitiva; ///esto sera un vector

	periodo=(double*)malloc(nt*sizeof(double));/// esto es un vector de n componentes. Es otra manera de crearlo
	primitiva=(double*)malloc(nt*sizeof(double));
	double cantidad[nt], Tmenos[nt], Tmas[nt];

	primitiva[0]=(b-a)*(funcion(a,u)+funcion(b,u))/2;
	periodo[0] = (4/w0)*primitiva[0];
	primitiva[1]=(b-a)*(funcion(a,u)+funcion(b,u)+2*funcion((a+(b-a)/2),u))/4;
	periodo[1] =(4/w0)*primitiva[1];
	primitiva[2]=(b-a)*(funcion(a,u)+funcion(b,u)+2*(funcion((a+(b-a)/3),u)+funcion((a+2*(b-a)/3),u)))/6 ;
	periodo[2] =(4/w0)*primitiva[2];

//Este bucle te calcula las primitivas y los periodos y los almacena en un array
	for (i=3; i<nt; i++){
            h_trapecio=(b-a)/(i+1);
			primitiva[i]=(h_trapecio/2)*(funcion(a,u)+funcion(b,u)+2*sumatorio(i,b,a,u));
			periodo[i]=(4/w0)*primitiva[i];
			for (j=1; j<i; j++){
				Tmenos[j]=fabs(periodo[j]-periodo[j-1]);
				Tmas[j]=periodo[j]+periodo[j+1];
				cantidad[j]=2*Tmenos[j]/Tmas[j];
				if (cantidad[j]<e){
					error=100*(periodo[j]-periodo_armonico_t)/periodo_armonico_t;
                    printf ("\nProceso detenido en la iteracion %d debido a que la tolerancia es:\n %.30lf\n", j+1, cantidad[j]);
					printf ("PRIMITIVA: %.15lf\n", primitiva[j]);
					printf ("PERIODO (s): %.15lf\n\n", periodo[j]);
					printf("PERIODO ARMONICO (s): %.15lf s\n", periodo_armonico_t);
					printf ("La variacion respecto el periodo armonico es: %lf %\n\n\n", error);
                    fprintf (archivo_trapecio, "\n\nProceso detenido en la iteración %d debido a que la tolerancia es:\n%.30lf\n", j+1, cantidad[j]);
					fprintf (archivo_trapecio, "\nPRIMITIVA: %.15lf\n", primitiva[j]);
					fprintf (archivo_trapecio, "PERIODO (s): %.15lf\n", periodo[j]);
					fprintf (archivo_trapecio, "PERIODO ARMÓNICO(s): %.15lf\n", periodo_armonico_t);
					fprintf (archivo_trapecio, "La variacion respecto el periodo armonico es: %lf\n\n\n", error);
					goto jump2;

				}
				else {
					if (j == nt-2){
						printf ("\n\nError. Tolerancia muy pequena.\nVuelve a iniciar el programa con mas pasos o un valor mayor para la tolerancia.\n\n");
						fprintf (archivo_trapecio, "\n\nError. Tolerancia muy pequena.\nVuelve a iniciar el programa con mas pasos o un valor mayor para la tolerancia.\n\n");
                    }
				}
			}
	}
	jump2:
    fclose (archivo_trapecio);
	printf ("\nEl proceso se ha concluido.\n\n");


};

double funcion(double fi,double u){
	double resultado;

	resultado = pow(sqrt(1-u*pow(sin(fi),2)),-1);

	return resultado;
}
double sumatorio(int nt,double b,double a,double u){
    double h_trapecio,suma;
    int j;

    h_trapecio=(b-a)/(nt+1);
    suma=0;
    for (j=1;j<nt+1;j++){
        suma=suma+funcion(a+j*h_trapecio,u);
        //printf ("%d: %lf\n", j, suma);
    }
    return suma;
}

void romberg(double fi0, double g, double L, double e, double u){
    FILE *archivo_romberg;
    int i, j, k,l, n, m;
	double sum, h, w0, periodo_armonico_r, error;
    archivo_romberg = fopen("calculo_romberg.txt", "w");

	printf ("Di hasta que orden N quieres calcular: ");
	scanf ("%d", &n);

	printf ("\n");

	a = 0;
	b = PI/2;
	h = b-a;
	w0 = sqrt(g/L);
	periodo_armonico_r=2*PI*sqrt(L/g);

    fprintf (archivo_romberg, "PROGRAMA REALIZADO POR:\n   Adrià Medeiros (NIU:1388403)\n   Luís Farfan (NIU:1388330)\n   Lucas Fuentealba(NIU:1388326)\n\n");
	fprintf (archivo_romberg, "CONDICIONES INICIALES:\n");
    fprintf (archivo_romberg, "Calculo hasta orden N: %d", n);
	fprintf (archivo_romberg, "\nCondicion inicial fi0: %lf", fi0);
	fprintf (archivo_romberg, "\nGravedad g: %lf", g);
	fprintf (archivo_romberg, "\nLongitud L: %lf", L);
	fprintf (archivo_romberg, "\nTolerancia e: %lf\n", e);

	double **R,*T;
	R=(double**)calloc(n,sizeof(double));
	for(i=0;i<=n;i++){
        R[i] = (double*)calloc(n,sizeof(double));
	}

	T=(double*)malloc(n*sizeof(double));

	double cantidad[n], Tmenos[n], Tmas[n];


	u=pow(sin(fi0/2),2);

	printf ("El proceso se esta realizando. Espera...\n\n");

	R[0][0]= 0.5*h*(funcion(0,u)+funcion(0.5*PI,u));// corresponde al metodo del trapecio base
	T[0]= 4/w0*R[0][0];

	for (i=1; i<n; i++){
		sum = 0;
		for (k=1; k<=(pow(2,i-1)); k++){//sumatorio que forma parte de la formula de la primera iteracion

			sum = sum + funcion((2*k-1)*h/pow(2,i), u);
			}
		R[i][0] = 0.5*(R[i-1][0] + (h/pow(2,i-1))*sum);

		for (j=1; j<=i;j++){
			R[i][j]= R[i][j-1]+(R[i][j-1]-R[i-1][j-1])/(pow(4,j)-1); //calcula la fila i, con i fijada
			if (i>=2){
                for(l=0;l<=j;l++){
                    T[l]=4/w0*R[i][l]; //calcula j-1 periodos de la fila i
                        for (m=1; m<l; m++){
                            Tmenos[m]=fabs(T[m]-T[m-1]);
                            Tmas[m]=T[m]+T[m+1];
                            cantidad[m]=2*Tmenos[m]/Tmas[m];
                            if (cantidad[m]<e){
                                error=100*(T[m]-periodo_armonico_r)/periodo_armonico_r;
                                printf ("\nProceso detenido en R[%d][%d] debido a que la tolerancia es:\n%.50lf\n", i, m, cantidad[m]);
                                printf ("El valor de la PRIMITIVA es: %.20lf\n", R[i][m]);
                                printf ("El valor del PERIODO es: %.20lf\n\n", T[m]);
                                printf("PERIODO ARMONICO (s): %.15lf s\n", periodo_armonico_r);
                                printf ("La variacion respecto el periodo armonico es: %lf %\n\n\n", error);
                                fprintf (archivo_romberg, "\nProceso detenido en R[%d][%d]", i, m);
                                fprintf (archivo_romberg, "\nTolerancia e =\n %.30lf\n\n", cantidad[m]);
                                fprintf (archivo_romberg, "\nPRIMITIVA: %.15lf\n", R[m][m]);
                                fprintf (archivo_romberg, "PERIODO (s): %.15lf\n", T[m]);
                                fprintf (archivo_romberg, "PERIODO ARMÓNICO(s): %.15lf\n", periodo_armonico_r);
                                fprintf (archivo_romberg, "La variacion respecto el periodo armonico es de: %lf %\n\n\n", error);
                                goto jump2;
                            }
                            else {
                                if (m==n-2){
                                    printf ("\n\nError. Tolerancia muy pequena.\nVuelve a iniciar el programa con mas pasos o un valor mayor para la tolerancia.\n\n");
                                fprintf (archivo_romberg, "\n\nError. Tolerancia muy pequena.\nVuelve a iniciar el programa con mas pasos o un valor mayor para la tolerancia.\n\n");
                                }
                            }
                        }
                }
            }
            else{
            }

        }
	}
    jump2:
	fclose (archivo_romberg);
	printf ("\nEl proceso se ha concluido.\n\n");


};
