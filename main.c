#include "Myheader.h"

double a, b, g, L, u, e, fi0;

int main()
{
    char myname[]="     Adria Medeiros (NIU:1388403)\n     Luis Farfan (NIU:1388330)\n     Lucas Fuentealba (NIU:1388326)\n";
    printf("Programa realizado por\n%s\n",myname);
    char repeticio;
    char procediment;


    do{
        printf ("\nDefine la condicion inicial fi0(en radianes): ");
        scanf ("%lf", &fi0);
        printf ("\nDa un valor para la gravedad g (m/s^2): ");
        scanf ("%lf", &g);
        printf ("\nDa un valor para la longitud L (m): ");
        scanf ("%lf", &L);
        printf ("\nDa un valor para la tolerancia: ");
        scanf ("%lf", &e);

        u = pow(sin(fi0/2),2);

        fflush(stdin);

        printf("\n\nElige que procedimiento quieres aplicar. \nPulsa 'T' para trapecio o 'R' para Romberg.\n");
        scanf("%c",&procediment);
        if ( procediment == 't' || procediment== 'T'){
            trapezi(fi0, g, L, e, u);
        }
        else{
            if (procediment== 'r' || procediment== 'R'){
                romberg(fi0, g, L, e, u);
            }
            else{
                printf("Por favor introduce un caracter valido\n");
            }
        }


        fflush(stdin);
        printf("Quieres volver a iniciar el programa?\n\nSi no quieres pulsa n \nSi  quieres pulsa cualquier otra tecla\n");
        scanf("%c",&repeticio);

    }while(repeticio !='n' && repeticio !='N');


    system("pause");
    return 0;
}
