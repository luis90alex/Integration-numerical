#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <locale.h>

#define PI 3.14159265358979323846

/// declarar todas las funciones

void trapezi(double fi0, double g, double L, double e, double u);
double funcion(double fi,double u);
double sumatorio(int n,double b,double a,double u);
double cant_periodo(int t);


void romberg(double fi0, double g, double L, double e, double u);
