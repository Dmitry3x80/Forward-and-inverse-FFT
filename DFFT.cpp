#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <random>

const double pi = 3.14159265358979323846;

// Function for performing fast forward discrete Fourier Transform (FFT)
std::vector<std::complex<double>> dft(const std::vector<std::complex<double>>& signal) {
    int N = signal.size();
    std::vector<std::complex<double>> spectrum(N);
    int T = 0; // Variable for storing the multiplicity value

    if (N % 5 == 0) {
        T = 5;
    }
    else if (N % 3 == 0) {
        T = 3;
    }
    else if (N % 2 == 0) {
        T = 2;
    }
    else {
        T = 1;
    }

    std::cout << "T:" << T << std::endl;

    for (int k = 0; k < N; ++k) {
        std::complex<double> sum(0.0, 0.0);
        for (int n = 0; n < T; ++n) {
            for (int m = 0; m < N / T; ++m) {
                double angle = - 2.0 * pi * k * (n + m * T) / N;
                std::complex<double> twiddle(std::cos(angle), std::sin(angle));
                sum += signal[n + m * T] * twiddle;
            }
        }
        spectrum[k] = sum;
    }

    return spectrum;
}

// Function for performing fast inverse discrete Fourier transform (IFFT)
std::vector<std::complex<double>> idft(const std::vector<std::complex<double>>& spectrum) {
    int N = spectrum.size();
    std::vector<std::complex<double>> signal(N);
    int T = 0; // Variable for storing the multiplicity value

    if (N % 5 == 0) {
        T = 5;
    }
    else if (N % 3 == 0) {
        T = 3;
    }
    else if (N % 2 == 0) {
        T = 2;
    }
    else {
        T = 1;
    }

    for (int k = 0; k < N; ++k) {
        std::complex<double> sum(0.0, 0.0);
        for (int n = 0; n < T; ++n) {
            for (int m = 0; m < N / T; ++m) {
                double angle = 2.0 * pi * k * (n + m * T) / N;
                std::complex<double> twiddle(std::cos(angle), std::sin(angle));
                sum += spectrum[n + m * T] * twiddle;
            }
        }
        signal[k] = sum * (1.0 / N);
    }

    return signal;
}

int main() {
    const int M = 10; // Size of an array of random complex numbers

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    std::vector<std::complex<double>> signal(M);

    for (int i = 0; i < M; ++i) {
        double real = dis(gen);
        double imag = dis(gen);
        signal[i] = std::complex<double>(real, imag);
    }

    std::vector<std::complex<double>> signal_copy(M);
    signal_copy = signal;

    // Output of the original array
    std::cout << "signal:" << std::endl;
    for (int i = 0; i < M; ++i) {
        std::cout << "signal[" << i << "]: " << signal[i] << std::endl;
    }

    int N = signal.size();

    std::vector<std::complex<double>> spectrum_dft = dft(signal);
    std::vector<std::complex<double>> spectrum_idft = idft(spectrum_dft);

    // Output of FFT results
    std::cout << "DFT spectrum:" << std::endl;
    for (int k = 0; k < N; ++k) {
        std::cout << "Spectrum[" << k << "]: " << spectrum_dft[k] << std::endl;
    }

    // Output of IFFT results
    std::cout << "IDFT spectrum:" << std::endl;
    for (int k = 0; k < N; ++k) {
        std::cout << "Spectrum[" << k << "]: " << spectrum_idft[k] << std::endl;
    }

    //Error Calculation
    std::vector<std::complex<double>> difference(signal_copy.size());
    std::vector<std::complex<double>> Error(signal_copy.size());

    //The difference between the original array and the array after transformations
    for (int i = 0; i < signal_copy.size(); ++i) {
        difference[i] = signal_copy[i] - spectrum_idft[i];
        Error[i] = difference[i] / signal_copy[i] * 100.0;
        std::cout << "difference[" << i << "]: " << difference[i] << std::endl;
    }

    //Percentage error
    for (int i = 0; i < signal_copy.size(); ++i) {
        difference[i] = signal_copy[i] - spectrum_idft[i];
        Error[i] = difference[i] / signal_copy[i] * 100.0;
        std::cout << "Error[" << i << "]: " << Error[i] << std::endl;
    }

    return 0;
}