#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include "Schrodinger.h"


using namespace std;

int main()
{
    int m, mt=0;//Número de veces que realizamos el experimento y veces que la particula ha cruzado el potencial
    int nd;//Tiempo tras el que aplicamos los detectores
    double K;//Coeficiente de transmisión
    double prob;//Variable auxiliar para almacenar la probabilidad
    double lambda,s,k0,norma;
    int ciclos,N;
    double *V,*modulo; //Potencial
    int j,l=0,t;//Indices el indice i queda reservado ya que nombro asi a la unidad imaginaria
    ofstream fich,fich2;//Fichero para almacenar los observables
    double mom,T,E;

    //Apertura de ficheros
    fich.open("Observables.txt");
    fich2.open("Energía.txt");

    //Número de veces que se va a repetir el experimento
    m=100;

    //Tiempo tras el que se detecta la particula
    nd=100;

    //Iniciacialización del generador de números aleatorios
    srand(time(NULL));

    //Introducción de N
    cout<<"Introduzca N"<<endl;
    cin>>N;

    //Vectores necesarios para el cálculo
    complex<double> phi[N];
    complex<double> a[N-1];
    complex<double> b[N-1];
    complex<double> chi[N];

    /** Recogida de parametros, inicialización de memoria dínamica y de algunas variables iniciales**/

    //Inicialización de la memoria dínamica

    V=new double [N];
    modulo=new double [N];

    //Introducción de la longitud de onda
    cout<<"Introduzca lambda"<<endl;
    cin>>lambda;



    //Introducción del número de ciclos
    do
    {
        cout<<"Introduzca el numero de ciclos a calcular entre 1 y "<<N/4 <<endl;
        cin>>ciclos;
    }
    while ((ciclos<1) || (ciclos >N/4));



    //Calculos necesarios para generar variables y condiciones
    k0=2.0*PI*ciclos/N;
    s=1.0/(4.0*k0*k0);

    //Repetición del experimento
    m=m+1;//Chapuza general para que la T quede comprendida entre 0 y 1
    genera: //Etiqueta para el goto
    while (l<m)
    {

        l++;

        //Generación de las funciones de onda iniciales y las alfas
        generador(s,V,k0,N,lambda,phi,a);

        //Calculo la norma inicial
        norma=0.0;
        for(j=0;j<N;j++)
            modulo[j]=norm(phi[j]);
        for(j=0;j<N;j++)
            norma+=modulo[j];

        //Normalizo la función de onda
        for(j=0;j<N;j++)
            phi[j]=phi[j]/sqrt(norma);

        t=0;
        /** Implementación del algoritmo para calcular la función de onda **/
        while (t<m*nd)
        {

            //Paso 2: Calculo de beta
            Beta(s,V,N,b,phi,a);

            //Paso 3: Calculo de chi
            Chi(N,a,b,chi);

            //Paso 4: Actualización de las funciones de onda
            Phi(N,chi,phi,modulo);

            t++;

            //Tras npasos hacemos las siguiente comparaciones
            if (t%nd==0)
            {



                prob=Pd(phi,N);

                if (aleatorio(prob)==false)
                {
                    mt++;
                    goto genera;

                }
                for (j=4*N/5;j<N;j++)
                    phi[j]=complex<double>(0.0,0.0);

                normalizador(phi,N);

                prob=Pi(phi,N);
                if (aleatorio(prob)==false)
                {
                    goto genera;
                }
                for (j=0;j<N/5;j++)
                    phi[j]=complex<double>(0.0,0.0);

                normalizador(phi,N);
            }
        }

    }

    mom=momento(phi,N);
    T=mom*mom;
    fich<<mom<<"\t"<<T<<endl;
    for (j=0;j<N;j++)
    {
        E=T+V[j];
        fich2<<j+1<<"\t"<<E<<endl;
    }

    K=1.0*mt/(1.0*m);

    cout<<"El valor del coeficiente de transmision es: "<<K<<endl;


    /** Borrado de la memoria dínamica y fin del programa **/


    //Borrado de memoria dínamica
    delete[] V;
    delete[] modulo;

    //Cierre de ficheros
    fich.close();


    return 0;
}
