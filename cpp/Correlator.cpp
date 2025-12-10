#include "Correlator.h"

#include <algorithm>
#include <numeric>

// Destructor
Correlator::~Correlator()
{
    if (plan_forward_a)
        fftw_destroy_plan(plan_forward_a);
    if (plan_forward_b)
        fftw_destroy_plan(plan_forward_b);
    if (plan_backward)
        fftw_destroy_plan(plan_backward);
}

// Calculate correlation
//! [in]  impulceResponce - Impulce Responce for Filters
//! [in]  signalIQ        - Singal to Correlate with
//! [out] corr_out        - Output correlation samples
//! return Output bits
std::vector<uint8_t> Correlator::correlate(const std::vector<std::vector<std::complex<double>>>& impulceResponce,
                                           const std::vector<std::complex<double>>&              signalIQ,
                                                 std::vector<std::vector<double>>&               corr_out)
{
    uint32_t numFilters = impulceResponce.size();
    uint32_t simbolSize = impulceResponce[0].size();
    uint32_t numSymbols = signalIQ.size() / simbolSize;

    std::vector<uint8_t> res(numSymbols * 2);

    // Find 4 correlations
    for (uint32_t filterId = 0; filterId < numFilters; ++filterId)
        correlateSignal(impulceResponce[filterId], signalIQ, corr_out[filterId]);

    uint32_t startSample = 0;
    std::vector<double> maxElements(numFilters);

    // Find symbols
    for (uint32_t iSymb = 0; iSymb < numSymbols; ++iSymb)
    {
        for (uint32_t filterId = 0; filterId < numFilters; ++filterId)
            maxElements[filterId] = *std::max_element(corr_out[filterId].begin() + startSample, corr_out[filterId].begin() + startSample + simbolSize);

        int maxIdx = std::max_element(maxElements.begin(), maxElements.end()) - maxElements.begin();

        // Decode sample
        switch (maxIdx)
        {
        case 0: // 00 simb
        {
            res[iSymb * 2 + 0] = 0;
            res[iSymb * 2 + 1] = 0;
            break;
        }
        case 1: // 01 simb
        {
            res[iSymb * 2 + 0] = 0;
            res[iSymb * 2 + 1] = 1;
            break;
        }
        case 2: // 10 simb
        {
            res[iSymb * 2 + 0] = 1;
            res[iSymb * 2 + 1] = 0;
            break;
        }
        case 3: // 11 simb
        {
            res[iSymb * 2 + 0] = 1;
            res[iSymb * 2 + 1] = 1;
            break;
        }
        default:
            std::cerr << "Error while decoding symbols: " << maxIdx << " id detected!" << std::endl;
            abort();
        }
        
        startSample += simbolSize;
    }

    return res;
}


// Process correlation with ouput data
//! [in]  data_a        - First signal to correlate
//! [in]  data_b        - Second signal to correlate
//! [out] corr_out      - Output correlation samples
void Correlator::correlateSignal(const std::vector<std::complex<double>>& data_a,
                                 const std::vector<std::complex<double>>& data_b,
                                       std::vector<double>&               corr_out)
{
    uint32_t size_a = data_a.size();
    uint32_t size_b = data_b.size();
    size_t size_out = size_a + size_b - 1;

    corr_out.resize(size_out);

    // std::complex<double> mean_a = std::accumulate(
    //     data_a.begin(), data_a.end(), std::complex<double>(0, 0)) / static_cast<double>(size_a);
    
    // std::complex<double> mean_b = std::accumulate(
    //     data_b.begin(), data_b.end(), std::complex<double>(0, 0)) / static_cast<double>(size_b);

    if (init)
    {
        a_centered.resize(size_a);
        b_centered.resize(size_b);
    }
    
    for (uint32_t i = 0; i < size_a; ++i)
        a_centered[i] = data_a[i];// - mean_a;
    
    for (uint32_t i = 0; i < size_b; ++i)
        b_centered[i] = data_b[i];// - mean_b;

    // double energy_a = 0.0;
    // double energy_b = 0.0;
    
    // for (const auto& val : a_centered)
    //     energy_a += std::norm(val);
    
    // for (const auto& val : b_centered)
    //     energy_b += std::norm(val);
    
    // double normalizer = std::sqrt(energy_a * energy_b);
    // if (normalizer < 1e-12)
    //     normalizer = 1e-12;

    if (init)
    {
        while (n_fft < size_out)
            n_fft <<= 1;
        a_fft.resize(n_fft);
        b_fft.resize(n_fft);
        corr_fft.resize(n_fft);
    }

    std::fill(a_fft.begin(), a_fft.end(), std::complex<double>{0,0});
    std::fill(b_fft.begin(), b_fft.end(), std::complex<double>{0,0});
    std::fill(corr_fft.begin(), corr_fft.end(), std::complex<double>{0,0});

    std::copy(a_centered.begin(), a_centered.end(), a_fft.begin());
    std::copy(b_centered.begin(), b_centered.end(), b_fft.begin());

    if (init)
    {
        plan_forward_a = fftw_plan_dft_1d(n_fft,
                                          reinterpret_cast<fftw_complex*>(a_fft.data()),
                                          reinterpret_cast<fftw_complex*>(a_fft.data()),
                                          FFTW_FORWARD,
                                          FFTW_ESTIMATE);
                                                
        plan_forward_b = fftw_plan_dft_1d(n_fft,
                                          reinterpret_cast<fftw_complex*>(b_fft.data()),
                                          reinterpret_cast<fftw_complex*>(b_fft.data()),
                                          FFTW_FORWARD,
                                          FFTW_ESTIMATE);
                                                
        plan_backward = fftw_plan_dft_1d(n_fft,
                                         reinterpret_cast<fftw_complex*>(corr_fft.data()),
                                         reinterpret_cast<fftw_complex*>(corr_fft.data()),
                                         FFTW_BACKWARD,
                                         FFTW_ESTIMATE);
        init = false;
    }

    fftw_execute(plan_forward_a);
    fftw_execute(plan_forward_b);

    for (uint32_t i = 0; i < n_fft; ++i)
        corr_fft[i] = a_fft[i] * std::conj(b_fft[i]);

    fftw_execute(plan_backward);

    for (uint32_t i = 0; i < size_out; ++i)
        corr_out[i] = std::abs(corr_fft[i]);//) / (n_fft * normalizer);

    return;
}
