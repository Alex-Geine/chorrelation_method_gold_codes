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
    std::vector<std::complex<double>> a_fft;
    std::vector<std::complex<double>> b_fft;
    std::vector<std::complex<double>> corr_fft;
    fftw_plan plan_forward_a;
    fftw_plan plan_forward_b;
    fftw_plan plan_backward;
    bool init = true;
    uint32_t n_fft = 1;
public: // functions

    // Calculate correlation
    //! [in]  impulceResponce - Impulce Responce for Filters
    //! [in]  signalIQ        - Singal to Correlate with
    //! [out] corr_out        - Output correlation samples
    //! return Output bits
    std::vector<uint8_t> correlate(const std::vector<std::vector<std::complex<double>>>& impulceResponce,
                                   const std::vector<std::complex<double>>&              signalIQ,
                                         std::vector<std::vector<double>>&               corr_out);




    // Destructor
    ~Correlator();

private: // functions

    // Process correlation with ouput data
    //! [in]  data_a        - First signal to correlate
    //! [in]  data_b        - Second signal to correlate
    //! [out] corr_out      - Output correlation samples
    void correlateSignal(const std::vector<std::complex<double>>& data_a,
                         const std::vector<std::complex<double>>& data_b,
                               std::vector<double>&               corr_out);
};

#endif //_CORRELATOR_H_
