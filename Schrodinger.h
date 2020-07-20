#ifndef SCHRODINGER_H_INCLUDED
#define SCHRODINGER_H_INCLUDED

#include <complex>
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>

using namespace std;

static const double PI= 3.14159265359; //Definicion de pi
static const double h=6.62607015e-34; //Constante de Planck
static const complex<double> i(0.0,1.0);//Definición de la unidad imaginaria
static const double hbarra=h/2*PI;//Definición de hbarra


void generador (double s, double V[], double k0,int N,double lambda,complex<double> phi[],complex<double>a[]);
void Beta (double s, double V[],int N,complex <double> b[],complex<double>phi[],complex<double>a[]);
void Chi (int N,complex<double>a[],complex<double>b[],complex<double>chi[]);
void Phi (int N,complex<double>chi[],complex<double>phi[],double modulo[]);

//Funciones necesarias para el voluntario
bool aleatorio (double Prob);
double Pd (complex<double> phi[], int N);
double Pi (complex<double> phi[], int N);
void normalizador (complex<double> phi[], int N);
double momento (complex<double> phi[], int N);


#endif // SCHRODINGER_H_INCLUDED
