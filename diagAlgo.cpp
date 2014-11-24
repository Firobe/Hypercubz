#include <iostream>
#include <cstdlib>
#include <vector>
#include <windows.h>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

long double coef( int , int , int );
vector<long double> eigenValues( int );
MatrixXd passage(int);
int binom(int, int);
void round(MatrixXd&, int);

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
    int l = pow(2, n);
    MatrixXd pas = passage(n);
    round(pas, l);
    MatrixXd invPas = pas.inverse();
    round(invPas, l);
    ///FILLING D WITH EIGENVALUES
    vector<long double> eigen;
    long double i = 0;
    while( 1 - ( 2.*i ) / ( double )n >= -1 ) {
        for(int k = 0 ; k < binom((int)i, n); k++)
            eigen.push_back( 1 - ( 2.*i ) / ( long double )n );
        i++;
    }
    MatrixXd d(l, l);
    for(int i = 0 ; i < l ; i++)
        for(int j = 0 ; j < l ; j++)
            if(i == j)
                d(i, j) = pow(eigen[i], p);
            else
                d(i, j) = 0;
    ///
    MatrixXd m = pas * d * invPas;
    round(m, l);
    cout << "Matrice diagonale :\n" << d << "\nMATRICE DE TRANSITION( " << pas.determinant() << " )\n" << pas << "\nINVERSE\n" << invPas << "\nRESULTAT\n" << m << endl;
    system( "pause" );
    return 0;
}
void round(MatrixXd& t, int l) {
    for(int i = 0 ; i < l ; i++)
        for(int j = 0 ; j < l ; j ++)
            if(t(i, j) < 0.00001 && t(i, j) > - 0.00001)
                t(i, j) = 0;
}

long double coef(const int i, const int j, const int n) {
    int a = 0, b = 0, l = pow(2, n);
    for(;;) {
        if( l <= 1)
            return 0;
        else if(i < a + l / 2 && j >= b + l / 2) {
            if(i - a == j - (b + l / 2))
                return 1. / (long double)n;
            else
                return 0.;
        } else if(i >= a + l / 2 && j < b + l / 2) {
            if(i - (a + l / 2) == j - b)
                return 1. / (long double)n;
            else
                return 0.;
        } else if(i < a + l / 2)
            l /= 2;
        else {
            l /= 2;
            a += l;
            b += l;
        }
    }
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

int binom(int k, int n) {
    int C[n + 1][k + 1];
    int i, j;
    for (i = 0; i <= n; i++) {
        for (j = 0; j <= min(i, k); j++) {
            if (j == 0 || j == i)
                C[i][j] = 1;
            else
                C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
        }
    }
    return C[n][k];
}

MatrixXd passage(int n) {
    int l = pow(2, n), c = 0;
    MatrixXd temp(l, l);
    vector<long double> spectre = eigenValues(n);
    cout << "Valeur propres : ";
    for(unsigned int i = 0 ; i < spectre.size() ; i++)
        cout << spectre[i] << " | ";
    cout << endl << "Calcul de la matrice de transition..." << endl;
    for(unsigned int v = 0 ; v < spectre.size() ; v++) {
        long double lambda = spectre[v];
        cout << endl << " -Valeur propre " << lambda;
        int m = binom(v, n);
        MatrixXd a(l - m, l - m);
        VectorXd b(l - m);
        for(int i = 0 ; i < l - m ; i++)
            for(int j = 0 ; j < l - m ; j++)
                if(i != j)
                    a(i, j) = coef(i, j, n);
                else
                    a(i, j) = coef(i, j, n) - lambda;
        for(int vec = 1 ; vec < m + 1 ; vec++) {
            for(int i = 0 ; i < l - m ; i++) {
                b(i) = 0;
                for(int k = l - 1 ; k > l - m - 1; k-- )
                    if(k == l - vec)
                        b(i) -= coef(i, k, n);
            }
            VectorXd x = a.partialPivLu().solve(b);
            cout << ".";
            for(int i = 0 ; i < l ; i ++)
                if(i < l - m)
                    temp(i, c) = x(i);
                else if(i == l - vec)
                    temp(i, c) = 1;
                else
                    temp(i, c) = 0;
            c++;
        }
    }
    cout << endl;
    return temp;
}
