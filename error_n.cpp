#include <iostream>
#include <fstream>
#include <complex>
#include <math.h>
#define MAX 200

using namespace std;

complex<double>* recur_FFT(complex<double>* arr, int N, int clockwise)
{
  if (N == 1) {
    return arr;
  }
  complex<double> *W;
  W = (complex<double> *)malloc(N * sizeof(complex<double>)); // 用複數存
  W[0] = 1; // W[0]: (1,0i) 初始值 x = 1, y = 0
  W[1] = polar(1., clockwise * 2. * M_PI / N); // mag 1, angle -2. * M_PI / N <- unit angle
  // initialize the rest of W
  for(int i = 2; i < N; i++)
    W[i] = pow(W[1], i);

  complex<double>* even = (complex<double>*)malloc(N / 2 * sizeof(complex<double>));
  complex<double>* odd = (complex<double>*)malloc(N / 2 * sizeof(complex<double>));
  for(int i = 0; i < N; i++) {
    if (i % 2 == 0) {
      even[i / 2] = arr[i];
    } else {
      odd[i / 2] = arr[i];
    }
  }

  complex<double> * even_result = recur_FFT(even, N/2, clockwise);
  complex<double> * odd_result = recur_FFT(odd, N/2, clockwise);

  complex<double> * result = (complex<double> *)malloc(N * sizeof(complex<double>));
  if (clockwise) {
  for(int i = 0; i < N / 2; i++) {
    result[i] = (even_result[i] + odd_result[i] * W[i]);
    result[i + N / 2] = (even_result[i] - odd_result[i] * W[i]);
  }
  } else {
    for(int i = N / 2 - 1; i >= 0; i--) {
      result[i] = (even_result[i] + odd_result[i] * W[i]);
      result[i - N / 2] = (even_result[i] - odd_result[i] * W[i]);
    }
  }

  return result;
};

double diff_square(complex<double> * original_coefficient, complex<double> * transformed_coefficient, int n, double x) {
    double sum = 0.0;
    for(int j = 0; j < n; j++) {
        double diff_real = (original_coefficient[j].real() - transformed_coefficient[j].real()) * pow(x, j);
        diff_real *= diff_real;
        double diff_imag = original_coefficient[j].imag() - transformed_coefficient[j].imag();
        diff_imag *= diff_imag;
        sum += diff_real;
        sum += diff_imag;
    }
    return sum;
}

int main()
{
    int n;
    double * x = new double[MAX];
    double * x_smaller_than1 = new double[MAX];
    int x_n;
    string x_data_number;

    complex<double> original_coefficient[MAX];
    string openfilename;
    
    cout << "Enter the amout of x value" << endl;
    cin >> x_n;

    

    return 0;
}
