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

/*
int main()
{
  int n;
  do {
    cout << "specify array dimension (MUST be power of 2)" << endl;
    cin >> n;
  } while(!check(n));
  double d;
  cout << "specify sampling step" << endl; //just write 1 in order to have the same results of matlab fft(.)
  cin >> d;
  complex<double> vec[MAX];
  cout << "specify the array" << endl;
  for(int i = 0; i < n; i++) {
    cout << "specify element number: " << i << endl;
    cin >> vec[i];
  }
  FFT(vec, n, d);
  cout << "...printing the FFT of the array specified" << endl;
  for(int j = 0; j < n; j++)
    cout << vec[j] << endl;
  return 0;
}
*/

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

    cout << "Enter the coefficient polynomial data of x" << endl;
    cin >> x_data_number;

    ifstream file("coefficient_polynomial" + x_data_number + ".txt");
    ifstream x_file("x_value.txt");
    for(int j = 0; j < x_n; j++) {
        x[j] = j;
    }
    for(int j = 0; j < x_n; j++) {
        x_smaller_than1[j] = pow(10,-1*j);
    }
    int surfix = 0;
    while(1)
    {
        
        file >> n;
        if( file.eof() ) break;
        string err_file_name = "error_data/error" + std::to_string(surfix);
        string err_small_1_file_name = "error_data/error_small" + std::to_string(surfix);
        ofstream error_file(err_file_name);
        ofstream error_small_file(err_small_1_file_name);

        for(int i = 0; i < n; ++i)
        {
            file >> original_coefficient[i];
        }

        for(int j = 0; j < n; j++)
          cout << original_coefficient[j] << endl;

        /* Sampling step, set to 1
        double d;
        cout << "specify sampling step" << endl; //just write 1 in order to have the same results of matlab fft(.)
        cin >> d;
        */

        complex<double> * ans = (complex<double> *)malloc(n * sizeof(complex<double>));
        ans = recur_FFT(original_coefficient, n, 1);
        // FFT(original_coefficient, n, 1);
        cout << "=================  FFT =============\n";
        for(int j = 0; j < n; j++)
        // cout << original_coefficient[j] << endl;
        cout << ans[j] << endl;

        cout << endl;
        cout << endl;
        cout << "=================  REVERSE =============\n";
        complex<double> * transformed_coefficient = (complex<double> *)malloc(n * sizeof(complex<double>));
        transformed_coefficient = recur_FFT(ans, n, -1);

        for(int j = 0; j < n; j++) {
        transformed_coefficient[j] = transformed_coefficient[j] / polar(n + 0.0, 0.0);
        }

        for(int j = 0; j < n; j++)
        cout << transformed_coefficient[j] << endl;

        cout << "================= ERROR RMS = sqrt(MSE) =============\n";
        double sum = 0.0;
        for(int j = 0; j < n; j++) {
            double diff_real = original_coefficient[j].real() - transformed_coefficient[j].real();
            diff_real *= diff_real;
            double diff_imag = original_coefficient[j].imag() - transformed_coefficient[j].imag();
            diff_imag *= diff_imag;
            sum += diff_imag;
        }
        double erro = sqrt(sum / 2 / n);
        cout << "error:" << erro << endl;

        cout << "================= ERROR X >= 1.0 =============\n";
        double error_sum = 0.0;
        double error_small_sum = 0.0;
        for(int i = 0; i < x_n; ++i)
        {
            error_sum = diff_square(original_coefficient, transformed_coefficient, n, x[i]);
            error_small_sum = diff_square(original_coefficient, transformed_coefficient, n, x_smaller_than1[i]);
            double error = sqrt(error_sum / n / 2);
            double error_small = sqrt(error_small_sum / n / 2);
            cout << "x = " << x[i] << "\t\t error = " << error << endl;
            cout << "x = " << x_smaller_than1[i] << "\t\t error = " << error_small << endl;
            error_file << x[i] << " " << error << endl;
            error_small_file << x_smaller_than1[i] << " " << error_small << endl;
        }

        error_file.close();
        error_small_file.close();
        surfix += 1;
    }

    return 0;
}
