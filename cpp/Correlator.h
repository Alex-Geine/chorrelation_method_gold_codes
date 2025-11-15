#ifndef _CORRELATOR_H_
#define _CORRELATOR_H_

#include <iostream>
#include <vector>
#include <complex>
#include <cstdint>
#include <fftw3.h>
#include <mutex>

class Correlator
{
private: // variables

    static std::mutex fftwMutex;  // Общий мьютекс для всех FFTW операций
    std::vector<std::complex<double>> a_centered;
    std::vector<std::complex<double>> b_centered;
    std::vector<std::complex<double>> a_fft;
    std::vector<std::complex<double>> b_fft;
    std::vector<std::complex<double>> corr_fft;
    fftw_plan plan_forward_a;
    fftw_plan plan_forward_b;
    fftw_plan plan_backward;
    bool init = true;
    uint32_t n_fft = 1;
private: // functions

// Calculate correlation
//! [in]  data_a        - First signal to correlate
//! [in]  data_b        - Second signal to correlate
//! [out] corr_out      - Output correlation samples
void findCorrelation(const std::vector<std::complex<double>>& data_a,
                            const std::vector<std::complex<double>>& data_b,
                                  std::vector<double>&               corr_out);


public: // functions

// Destructor
~Correlator();

// Process correlation with ouput data
//! [in]  data_a        - First signal to correlate
//! [in]  data_b        - Second signal to correlate
//! [out] corr_out      - Output correlation samples
//! [out] max_metric_id - Index of maximum metric of the correlation
void correlate(const std::vector<std::complex<double>>& data_a,
                      const std::vector<std::complex<double>>& data_b,
                            std::vector<double>&               corr_out,
                            uint32_t&                          max_metric_id);

// Process correlation without ouput data
//! [in]  data_a        - Output generated AWGN
//! [in]  data_b        - Size of the output data
//! [out] max_metric_id - Index of maximum metric of the correlation
void correlate(const std::vector<std::complex<double>>& data_a,
                      const std::vector<std::complex<double>>& data_b,
                            uint32_t&                          max_metric_id);

};

#endif //_CORRELATOR_H_
