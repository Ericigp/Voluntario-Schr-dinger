/*********************************************************/
/**  Archivo que contiene las funciones necesarias para **/
/**  implementar el algoritmo de resolución de la       **/
/**  ecuación de Schrödinger.                           **/
/**  Las funciones se llaman igual que las variables    **/
/**  que calculan solo que con mayuscula.               **/
/**  Esto lo hago pq tengo 0 originalidad para poner    **/
/**  nombres.                                           **/
/*********************************************************/

#include "Schrodinger.h"

using namespace std;

void generador (double s, double V[], double k0,int N,double lambda,complex<double> phi[],complex<double>a[])
{
    int j;
    complex<double>Phi;

    //Condiciones de contorno
    phi[0]=phi[N-1]=complex<double>(0.0,0.0);
    V[0]=V[N-1]=0.0;
    a[N-1]=complex<double>(0.0);

    //Calculo de la función de onda
    for(j=1;j<(N-1);j++)
    {
        phi[j]=exp(j*k0*i)*exp(-8.0*(4.0*j-N)*(4.0*j-N)/(N*N*1.0));
    }

    //Inicialización del potencial
    for (j=0;j<N;j++)
    {
        if ((j>=(2.0*N/5.0)) && (j<=(3.0*N/5.0)))
            V[j]=(lambda*k0*k0);
        else
            V[j]=0.0;
    }

    //Inicialización de alfa
    for(j=(N-2);j>0;j--)
    {
        a[j]=-1.0/(-2.0-V[j+1]+a[j+1]+(2.0*i)/s);
    }

    return;
}

void Beta (double s, double V[],int N,complex <double> b[],complex<double>phi[],complex<double>a[])
{
    int j;

    //Condiciones de contorno
    b[N-1]=complex<double>(0.0,0.0);

    for(j=(N-2);j>0;j--)
    {
        b[j]=(((4.0*i*phi[j+1])/s)-b[j+1])/(-2.0-V[j+1]+a[j+1]+(2.0*i)/s);
    }

}

void Chi (int N,complex<double>a[],complex<double>b[],complex<double>chi[])
{
    int j;

    //Condiciones de contorno
    chi[0]=complex<double>(0.0,0.0);
    chi[N-1]=complex<double>(0.0,0.0);

    for(j=0;j<(N-2);j++)
    {
        chi[j+1]=a[j]*chi[j]+b[j];
    }
    return;
}

void Phi (int N,complex<double>chi[],complex<double>phi[],double modulo[])
{
    int j;
    for(j=0;j<N;j++)
    {
        phi[j]=chi[j]-phi[j];
        modulo[j]=norm(phi[j]);
    }
    return;
}

//Para hacer una representación del potencial dinamico he de hacer esta chapuza
void Script(string pot2)
{
    ifstream entrada;
    ofstream salida;
    string texto [100];
    int i=0,n;

    entrada.open("Script probabilidadog.txt");
    salida.open("Script probabilidad.txt");

    while (!entrada.eof())
    {
        getline(entrada,texto[i]);
        i++;
    }

    n=i;
    texto[13]="pot="+pot2;

    for (i=0;i<n;i++)
        salida<<texto[i]<<endl;

    entrada.close();
    salida.close();

    return;
}

bool aleatorio (double Prob)
{
    double x;
    x=1.0*rand()/RAND_MAX;

    if (x<Prob)
        return false;
    else
        return true;
}

double Pd (complex<double> phi[], int N)
{
    double suma=0.0;
    int i;

    for (i=(4*N/5);i<N;i++)
        suma+=norm(phi[i]);

    return suma;
}

double Pi (complex<double> phi[], int N)
{
    double suma=0.0;
    int i;

    for (i=0;i<(N/5);i++)
        suma+=norm(phi[i]);

    return suma;
}

void normalizador (complex<double> phi[], int N)
{
    double suma=0.0;
    int i;

    for(i=0;i<N;i++)
        suma+=norm(phi[i]);

    suma=sqrt(suma);

    for (i=0;i<N;i++)
        phi[i]/=suma;
}

double momento (complex<double> phi[], int N)
{
    complex<double> suma(0.0,0.0);
    int j;

    for(j=0;j<N-1;j++)
    {
        suma+=conj(phi[j])*(phi[j+1]-phi[j-1])*0.5;
    }

    return real(-i*suma);

}
