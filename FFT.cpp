#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

#define ld long double

using namespace std;

struct complex{
    ld real;
    ld imaginary;
    complex(){ real = 0, imaginary = 0; }
    complex(ld r, ld i){ real = r, imaginary = i; }
    void setWithPolarVal(ld r, ld theta){
        real = r * cos(theta);
        imaginary = r * sin(theta);
    }
    complex operator*(const complex& a) const{
        return complex(
            this->real * a.real - this->imaginary * a.imaginary,
            this->real * a.imaginary + this->imaginary * a.real
        );
    }
    complex operator+(const complex& a) const{
        return complex(
            this->real + a.real,
            this->imaginary + a.imaginary
        );
    }
    complex operator-(const complex& a) const{
        return complex(
            this->real - a.real,
            this->imaginary - a.imaginary
        );
    }
    complex operator/(const int& a) const{
        return complex(
            this->real / a,
            this->imaginary / a
        );
    }
};
vector<complex> fftRec(vector<complex>& pol, complex w){
    int n = pol.size();
    if(n == 1) return pol;
    vector<complex> APol;
    vector<complex> BPol;
    for(int i = 0; i < n; i++){
        if(i % 2) BPol.push_back(pol[i]);
        else APol.push_back(pol[i]);
    }
    vector<complex> A = fftRec(APol, w * w);
    vector<complex> B = fftRec(BPol, w * w);
    vector<complex> res(n);
    complex ws(1,0);
    int k = n/2;
    for(int i = 0; i < k; i++){
        res[i] = A[i] + ws * B[i];
        res[i+k] = A[i] - ws * B[i];
        ws = ws * w;
    }
    return res;
}
int nextPow2(int n){
    // if N is a power of two simply return it
    if (!(n & (n - 1))) return n;
    // else set only the left bit of most significant bit
    return 0x8000000000000000 >> (__builtin_clzll(n) - 1);
}
vector<complex> setupExpression(vector<complex>& pol){
    int n = pol.size();
    vector<complex> res(nextPow2(n));
    for(int i = 0; i < n; i++) res[i] = pol[i];
    return res;
}
vector<complex> FFT(vector<complex>& pol){
    pol = setupExpression(pol);
    int n = pol.size();
    complex w;
    w.setWithPolarVal(1, 2 * M_PI / n);
    return fftRec(pol, w);
}
vector<complex> IFFT(vector<complex>& val){
    val = setupExpression(val);
    int n = val.size();
    complex w;
    w.setWithPolarVal(1, -2 * M_PI / n);
    vector<complex> res = fftRec(val, w);
    for(int i = 0; i < n; i++){
        res[i] = res[i] / n;
    }
    return res;
} 
int main(){

    vector<complex> pol = {
        complex(2,0),
        complex(1,0),
        complex(-1,0),
        complex(5,0),
        complex(0,0),
        complex(3,0),
        complex(0,0),
        complex(-4,0),
        complex(10,0)
    };
    complex w(0,1);
    int n = pol.size();

    vector<complex> _fft = FFT(pol);
    for(int i = 0; i < n; i++){
        cout << _fft[i].real << " " << _fft[i].imaginary << endl;
    }

    cout << "---" << endl;

    vector<complex> _ifft = IFFT(_fft);
    for(int i = 0; i < n; i++){
        cout << fixed << setprecision(2) << _ifft[i].real << endl;
    }

    return 0;
}