#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <complex>

#define ld long double
#define complex complex<double>

using namespace std;

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
    complex w = polar(1.0, 2 * M_PI / n);
    return fftRec(pol, w);
}
vector<complex> IFFT(vector<complex>& val){
    val = setupExpression(val);
    int n = val.size();
    complex w = polar(1.0, -2 * M_PI / n);
    vector<complex> res = fftRec(val, w);
    for(int i = 0; i < n; i++){
        res[i] /= complex(n,0);
    }
    return res;
} 

vector<vector<complex>> Transpone(vector<vector<complex>>& m){
    vector<vector<complex>> res(m[0].size(), vector<complex>(m.size()));
    for(int i = 0; i < m.size(); i++){
        for(int j = 0; j < m[i].size(); j++){
            res[j][i] = m[i][j];
        }
    }
    return res;
}
vector<vector<complex>> Setup2D(vector<vector<complex>>& m){
    int n = nextPow2(m.size());
    int l = nextPow2(m[0].size());
    vector<vector<complex>> res(n, vector<complex>(l));
    for(int i = 0; i < m.size(); i++){
        for(int j = 0; j < m[i].size(); j++){
            res[i][j] = m[i][j];
        }
    }
    return res;
}
vector<vector<complex>> ApplyFFTToRows(vector<vector<complex>>& m, complex w){
    for(int i = 0; i < m.size(); i++){
        complex tmp = w;
        m[i] = fftRec(m[i], tmp);
    }
    return m;
}
vector<vector<complex>> FFT2D(vector<vector<complex>>& m){
    m = Setup2D(m);

    int n = m[0].size();
    complex w = polar(1.0, 2 * M_PI / n);
    m = ApplyFFTToRows(m,w);
    m = Transpone(m);

    n = m[0].size();
    w = polar(1.0, 2 * M_PI / n);
    m = ApplyFFTToRows(m,w);
    m = Transpone(m);

    return m;
}
vector<vector<complex>> IFFT2D(vector<vector<complex>>& m){
    m = Setup2D(m);

    int n = m[0].size();
    complex w = polar(1.0, -2 * M_PI / n);
    m = ApplyFFTToRows(m,w);
    m = Transpone(m);
    for(int i = 0; i < m.size(); i++){
        for(int j = 0; j < m[0].size(); j++){
            m[i][j] /= complex(n,0);
        }
    }

    n = m[0].size();
    w = polar(1.0, -2 * M_PI / n);
    m = ApplyFFTToRows(m,w);
    m = Transpone(m);
    for(int i = 0; i < m.size(); i++){
        for(int j = 0; j < m[0].size(); j++){
            m[i][j] /= complex(n,0);
        }
    }

    return m;
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
        cout << _fft[i] << endl;
    }

    cout << "---" << endl;

    vector<complex> _ifft = IFFT(_fft);
    for(int i = 0; i < n; i++){
        cout << fixed << setprecision(2) << _ifft[i] << endl;
    }

    cout << endl << endl << endl;

    vector<vector<complex>> m = {
        {complex(1.0,0), complex(2.0,0), complex(-3.0,0), complex(5.0,0)},
        {complex(-17.0,0), complex(20.0,0), complex(0.0,0), complex(0.0,0)},
        {complex(0.0,0), complex(1.0,0), complex(-3.0,0), complex(23.0,0)},
        {complex(6.0,0), complex(-7.0,0), complex(8.0,0), complex(-9.0,0)}
    };

    m = FFT2D(m);

    for(int i = 0; i < m.size(); i++){
        for(int j = 0; j < m[i].size(); j++){
            cout << m[i][j] << " ";
        }
        cout << endl;
    }

    cout << "---" << endl;

    m = IFFT2D(m);

    for(int i = 0; i < m.size(); i++){
        for(int j = 0; j < m[i].size(); j++){
            cout << m[i][j] << " ";
        }
        cout << endl;
    }

    return 0;
}