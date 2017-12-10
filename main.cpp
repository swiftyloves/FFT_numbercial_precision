#include <iostream>
#include <complex>
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

int main()
{
    int n;
    complex<double> vec[MAX];

    ofstream myfile;
    myfile.open ("mse_double.csv");

    ifstream file("coefficient.txt");
    while(1)
    {
        file >> n;
        if( file.eof() ) break;

        for(int i = 0; i < n; ++i)
        {
            file >> vec[i];
        }

        for(int j = 0; j < n; j++)
          cout << vec[j] << endl;

        /* Sampling step, set to 1
        double d;
        cout << "specify sampling step" << endl; //just write 1 in order to have the same results of matlab fft(.)
        cin >> d;
        */

        complex<double> * ans = (complex<double> *)malloc(n * sizeof(complex<double>));
        ans = recur_FFT(vec, n, 1);
        // FFT(vec, n, 1);
        cout << "=================  FFT =============\n";
        for(int j = 0; j < n; j++)
        // cout << vec[j] << endl;
        cout << ans[j] << endl;

        cout << endl;
        cout << endl;
        cout << "=================  REVERSE =============\n";
        complex<double> * reverse_ans = (complex<double> *)malloc(n * sizeof(complex<double>));
        reverse_ans = recur_FFT(ans, n, -1);

        for(int j = 0; j < n; j++) {
        reverse_ans[j] = reverse_ans[j] / polar(n + 0.0, 0.0);
        }

        for(int j = 0; j < n; j++)
        cout << reverse_ans[j] << endl;
    }


    // myfile << "Writing this to a file.\n";
    myfile.close();

    return 0;
}
