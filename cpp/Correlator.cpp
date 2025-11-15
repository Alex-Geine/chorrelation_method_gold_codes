#include "Correlator.h"

#include <algorithm>
#include <numeric>

std::mutex Correlator::fftwMutex;

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
//! [in]  data_a        - First signal to correlate
//! [in]  data_b        - Second signal to correlate
//! [out] corr_out      - Output correlation samples
void Correlator::findCorrelation(const std::vector<std::complex<double>>& data_a,
                                 const std::vector<std::complex<double>>& data_b,
                                 std::vector<double>&                     corr_out)
{
    uint32_t size_a = data_a.size();
    uint32_t size_b = data_b.size();
    size_t size_out = size_a + size_b - 1;

    if (init)
        corr_out.resize(size_out);

    std::complex<double> mean_a = std::accumulate(
        data_a.begin(), data_a.end(), std::complex<double>(0, 0)) / static_cast<double>(size_a);
    
    std::complex<double> mean_b = std::accumulate(
        data_b.begin(), data_b.end(), std::complex<double>(0, 0)) / static_cast<double>(size_b);

    if (init)
    {
        a_centered.resize(size_a);
        b_centered.resize(size_b);
    }
    
    for (uint32_t i = 0; i < size_a; ++i)
        a_centered[i] = data_a[i] - mean_a;
    
    for (uint32_t i = 0; i < size_b; ++i)
        b_centered[i] = data_b[i] - mean_b;

    double energy_a = 0.0;
    double energy_b = 0.0;
    
    for (const auto& val : a_centered)
        energy_a += std::norm(val);
    
    for (const auto& val : b_centered)
        energy_b += std::norm(val);
    
    double normalizer = std::sqrt(energy_a * energy_b);
    if (normalizer < 1e-12)
        normalizer = 1e-12;

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

    std::lock_guard<std::mutex> lock(fftwMutex);

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
        corr_out[i] = std::abs(corr_fft[i]) / (n_fft * normalizer);

    return;
}

// Process correlation with ouput data
//! [in]  data_a        - First signal to correlate
//! [in]  data_b        - Second signal to correlate
//! [out] corr_out      - Output correlation samples
//! [out] max_metric_id - Index of maximum metric of the correlation
void Correlator::correlate(const std::vector<std::complex<double>>& data_a,
                           const std::vector<std::complex<double>>& data_b,
                                 std::vector<double>&               corr_out,
                                 uint32_t&                          max_metric_id)
{
    if (data_b.size() > data_a.size())
        throw std::runtime_error("Error in correlate function. Size of data_a less then the size of data_b");

    findCorrelation(data_a, data_b, corr_out);

    // Находим итератор на максимальный элемент
    int peak_index = std::max_element(corr_out.begin(), corr_out.end()) - corr_out.begin();
    
    std::cout << "max el: " << corr_out[peak_index] << std::endl;
    // Вычисляем индекс
    max_metric_id = peak_index;// - (data_b.size() - 1);

    return;
}

// Process correlation with ouput data
//! [in]  data_a        - First signal to correlate
//! [in]  data_b        - Second signal to correlate
//! [out] max_metric_id - Index of maximum metric of the correlation
void Correlator::correlate(const std::vector<std::complex<double>>& data_a,
                           const std::vector<std::complex<double>>& data_b,
                                 uint32_t&                          max_metric_id)
{
    if (data_b.size() > data_a.size())
        throw std::runtime_error("Error in correlate function. Size of data_a less then the size of data_b");

    std::vector<double> corr_out;

    findCorrelation(data_a, data_b, corr_out);

    // Находим итератор на максимальный элемент
    int peak_index = std::max_element(corr_out.begin(), corr_out.end()) - corr_out.begin();

    
    // Вычисляем индекс
    max_metric_id = peak_index;// - (data_b.size() - 1);

    return;
}