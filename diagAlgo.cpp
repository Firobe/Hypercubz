#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <windows.h>
#include <cmath>
#include <chrono>

using namespace std;

vector<long double> eigenValues( int );
long double coefPas(int i, int j, int n);
int sum(int i, int n);

int main() {
    int n, p;
    cout << "Calcul optimise des puissances d'une matrice de transition d'hypercube" << endl
         <<  "------------------------------------------------------------------------" << endl;
    do {
        cout <<  "Dimension de l'hypercube : ";
        cin >> n;
    } while( n < 0 );
    do {
        cout << "Nombre de puissances a effectuer : ";
        cin >> p;
    } while( p < 1 );
    auto start = chrono::high_resolution_clock::now();
    int l = pow(2, n);
    vector<long double> eigen = eigenValues(n);
    double temp;
    for(int j = 0 ; j < l ; j++) {
        temp=0;
        for(int k = 0 ; k < l ; k++)
            temp += pow(eigen[sum(k, n)], p) * coefPas(0, k, n) * coefPas(j, k, n);
            if(temp < 0.00001 && temp > - 0.00001)
                temp = 0;
        cout << temp << " ";
        }
    auto finish = chrono::high_resolution_clock::now();
    cout << "Operation effectuée en " << (double)chrono::duration_cast<std::chrono::microseconds>(finish-start).count()/1000000 << " s." << endl;
    system("pause");
    return 0;
}

vector<long double> eigenValues( int n ) {
    vector<long double> temp;
    long double i = 0;
    while( 1 - ( 2.*i ) / (long double )n >= -1 ) {
        temp.push_back( 1 - ( 2.*i ) / (long  double )n );
        i++;
    }
    return temp;
}

long double coefPas(int i, int j, int n) {
    int scal = 0;
    for(int k = 0 ; k < n ; k++) {
        scal += (i % 2) * (j % 2);
        i = (int)floor(i / 2);
        j = (int)floor(j / 2);
    }
    long double p = -(long double)n / 2.;
    return pow(-1, scal) * pow(2, p);
}

int sum(int i, int n)
{
    int sum = 0;
    for(int k = 0 ; k < n ; k++) {
        sum+= i%2;
        i = (int)floor(i/2);
    }
    return sum;
}
