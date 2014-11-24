#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>

using namespace std;
bool conv(const vector<double>&, int);
int opti(const vector<double>&);
double getPas(const int, const int, const int);
vector<double> product (vector<double>&, vector<double>&);
int zeros = 0;

template <typename T>
vector<T> operator+ (const vector<T> &lhs, const vector<T> &rhs)
{
	vector<T> ret (lhs);
	ret.insert (ret.end(), rhs.begin(), rhs.end());
	return ret;
}
template <typename T> struct Plus
{
	vector<T> operator() (const vector<T> &lhs, const vector<T> &rhs){return lhs+rhs;}
};
template <typename T> void print (const vector<T> &v)
{
	for(const T &i : v)
		cout << i << ' ';
}
vector<double> add(const vector<double>& a, const vector<double>& b)
{
	vector<double> ret(a.size());
	for(unsigned int i = 0 ; i < a.size() ; i++)
		ret[i]=a[i]+ b[i];
	return ret;
}
int main()
{
	int n, pui=0;
	bool co;
	cout << "Recherche de convergence ? ";
	cin >> co;
	cout << "Dimension ? ";
	cin >> n;
	if(!co){
		cout << "Nb de puissances ? ";
		cin >> pui;}
	auto start = chrono::high_resolution_clock::now();
	vector<double> passage (pow(2,n));
	vector<double> calc (pow(2,n));
	for(int i = 0 ; i<pow(2,n) ; i++)
		calc[i]=passage[i]=getPas(0,i,n);
	if(!co){
		for(int i = 0; i<pui-1 ;i++)
			calc = product(calc, passage);
		print(calc);}
	else{
		while(!conv(calc, n)){
			pui++;
			calc = product(calc, passage);}
		cout << "La matrice a convergé après "<<pui<<" puissances."<<endl;}
	auto finish = chrono::high_resolution_clock::now();
	cout << "Operation effectuée en " << (double)chrono::duration_cast<std::chrono::microseconds>(finish-start).count()/1000000 << " s." << endl;
	cout << "Optimisations : " << zeros << endl;
	return 0;
}

vector<double> product (vector<double>& a, vector<double>& b)
{
	double piv = b[b.size()/2];
	if(a.size()==1)
		return (vector<double>){a[0]*b[0]};
	int optiA = opti(a), optiB = opti(b);
	if(optiA == 0 || optiB == 0){
		vector<double> ret(a.size(),0);
		return ret;}
	if(optiA == 1){
		zeros++;
		vector<double> ret(a.size());
		for(int i=0;i<ret.size();i++)
			ret[i]=b[i]*a[0];
		return ret;}
	if(optiB == 1){
		zeros++;
		vector<double> ret(a.size());
		for(int i=0;i<ret.size();i++)
			ret[i]=a[i]*b[0];
		return ret;}
	vector<double> a1 ( a.begin(), a.begin()+(a.size()/2));
	vector<double> a2 ( a.begin()+(a.size()/2), a.end());
	vector<double> b1 ( b.begin(), b.begin()+(b.size()/2));
	vector<double> c1(a1.size());
	vector<double> c2(a1.size());
	for(int i = 0; i<c1.size();i++){
		c1[i] = a2[i]*piv;
		c2[i] = a1[i]*piv;}
	return add(product(a1,b1), c1) + add(c2, product(a2, b1));
}

double getPas(const int i, const int j, const int n)
{
	int a = 0, b = 0, l = pow(2,n);
	for(;;) {
		if( l<=1)
			return 0;
		else if(i<a+l/2 && j>=b+l/2){
			if(i-a==j-(b+l/2))
				return 1./(float)n;
			else
				return 0;}
		else if(i>=a+l/2 && j<b+l/2){
			if(i-(a+l/2)==j-b)
				return 1./(float)n;
			else
				return 0;}
		else if(i<a+l/2)
			l/=2;
		else{
			l/=2;
			a+=l;
			b+=l;}
	}
}

bool conv(const vector<double>& v, const int n)
{
	double sum=0;
	for(unsigned int i=0;i<v.size();i++)
		if(v[i]!=0 && !(v[i]<=1./(pow(2,n-1))+0.0000000001 && v[i]>=1./(pow(2, n-1))-0.0000000001))
			return false;
	return true;
}

int opti(const vector<double>& v)
{
	for(int i=1;i<v.size();i++)
		if(v[i]!=0)
			return -1;
	if(v[0]==0)
		return 0;
	return  1;
}

