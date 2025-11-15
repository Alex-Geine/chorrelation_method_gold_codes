#include "Generator.h"
#include "Correlator.h"
#include <fstream>
#include <string>

constexpr double PI_2 = 6.28318530718;

namespace Utils
{
template <typename T>
void print(std::vector<T> data, std::string name)
{
    std::cout << name << ", size: " << data.size() << std::endl;
    for (auto& it: data)
        std::cout << it << " ";
    std::cout << std::endl;
}

template <typename T>
void write(std::vector<T> data, std::string name)
{
    std::ofstream outFile;
    outFile.open(name);

    if (!outFile.is_open())
        throw std::runtime_error("Cannot open file: " + name);

    for (uint32_t i = 0; i < data.size() - 1; ++i)
       outFile << data[i] << ", ";

    outFile << data[data.size() - 1] << std::endl;

    outFile.close();
}

void writeBer(double val, std::string name)
{
    std::ofstream outFile;
    outFile.open(name, std::ios::app);

    if (!outFile.is_open())
        throw std::runtime_error("Cannot open file: " + name);

    outFile << val << "\n";

    outFile.close();
}

}; // Utils

// Configure function
//! [in] std  - Standart noise deviation
//! [in] mean - Mean value 
void RandomGenerator::config(double std, double mean)
{
    dist = std::normal_distribution<double>(mean, std);
}

// Generate AWGN
//! [out] data_out  - Output generated AWGN
//! [in]  size      - Size of the output data
void RandomGenerator::generateAwgn(std::vector<std::complex<double>>& data_out, uint32_t size)
{
    if (size == 0)
        throw std::runtime_error("Error in generateAwgn function! Unsupported data size: " + size);

    data_out.resize(size);

    for (uint32_t i = 0; i < size; ++i)
        data_out[i] = std::complex<double>(dist(gen), dist(gen));

    return;
}

// Function for generating shifted in TD signal
//! [in]  sample_freq      - Sample frequency of the signal
//! [in]  d_t              - Time offset in sec
//! [in]  shifted_size     - Size of shifted signal
//! [in]  data_in          - Input data samples
//! [out] data_out         - Output shifted signal
void SignalGenerator::generateShiftedSignal(double sample_freq, double d_t, uint32_t shifted_size,
                                            const std::vector<std::complex<double>>& data_in,
                                                  std::vector<std::complex<double>>& data_out)
{
    if (sample_freq <= 0)
        throw std::runtime_error("Error in SignalGenerator::generateShiftedSignal."
                                 " Invalid sample_freq: " + std::to_string(sample_freq));
    if (d_t < 0)
        throw std::runtime_error("Error in SignalGenerator::generateShiftedSignal."
                                 " Invalid d_t: " + std::to_string(d_t));

    uint32_t size = data_in.size();

    if (size < shifted_size)
        throw std::runtime_error("Error in SignalGenerator::generateShiftedSignal."
                                 " Invalid shifted_size: " + std::to_string(shifted_size) +
                                 std::string("while size of data_in: ") + std::to_string(size));

    // Calculating n_shift  
    uint32_t n_shift = d_t;// * sample_freq;

    std::cout << "n_shift: " << n_shift << std::endl;

    if (n_shift + shifted_size > size)
        throw std::runtime_error("Error in SignalGenerator::generateShiftedSignal." +
                                 std::string(" Invalid n_shift: ") + std::to_string(n_shift) +
                                 std::string(", n_shift is size / fd / d_t, where size: ") +
                                 std::to_string(size) + std::string(", fd: ") + std::to_string(sample_freq) +
                                 std::string(", d_t: ") + std::to_string(d_t));

    if (data_out.empty())
        data_out.resize(shifted_size);

    auto iter_data_in = data_in.begin() + n_shift;

    std::copy(iter_data_in, iter_data_in + shifted_size, data_out.begin());

    return;
}

// Function for generating shifted in TD signal
double SignalGenerator::generateShiftedSignal(double sample_freq, uint32_t shifted_size,
                             const std::vector<std::complex<double>>& data_in,
                             double seed,
                             std::vector<std::complex<double>>& data_out)
{
    double dt = 0;

    if (sample_freq <= 0)
        throw std::runtime_error("Error in SignalGenerator::generateShiftedSignal."
                                 " Invalid sample_freq: " + std::to_string(sample_freq));

    uint32_t size = data_in.size();

    if (size < shifted_size)
        throw std::runtime_error("Error in SignalGenerator::generateShiftedSignal."
                                 " Invalid shifted_size: " + std::to_string(shifted_size) +
                                 std::string("while size of data_in: ") + std::to_string(size));

    dt = (size - shifted_size) * seed;

    // Calculating n_shift  
    uint32_t n_shift = dt;// * sample_freq;

    std::cout << "n_shift: " << n_shift << std::endl;

    if (n_shift + shifted_size > size)
        throw std::runtime_error("Error in SignalGenerator::generateShiftedSignal." +
                                 std::string(" Invalid n_shift: ") + std::to_string(n_shift) +
                                 std::string(", n_shift is size / fd / d_t, where size: ") +
                                 std::to_string(size) + std::string(", fd: ") + std::to_string(sample_freq) +
                                 std::string(", d_t: ") + std::to_string(dt));

    if (data_out.empty())
        data_out.resize(shifted_size);

    auto iter_data_in = data_in.begin() + n_shift;

    std::copy(iter_data_in, iter_data_in + shifted_size, data_out.begin());

    return dt;
}


// Add noise in data
//! [in/out] data     - Input/Output data
//! [in]     snr      - Signal to Noise Ratio
void NoiseInjector::addNoise(std::vector<std::complex<double>>& data, double snr)
{
    double energy = 0;

    for (auto& it: data)
        energy += std::abs(it * std::conj(it));

    // snr = log(energy / noise)
    // 10^snr = enegry / noise
    // noise = energy / 10^snr
    double lin_snr = std::pow(10, snr);
    double noise_power = energy / lin_snr;

    // need to process mean energy
    // Very big questions
    gen.config(std::sqrt(noise_power / 2.));

    uint32_t size = data.size();

    std::vector<std::complex<double>> noise;
    gen.generateAwgn(noise, size);

    for (uint32_t i = 0; i < size; ++i)
        data[i] += noise[i];

    return;
}

//! Configurate signal generator
//! [in] params - Configuration parameters
void BaseGenerator::configure(const cfg& params)
{
    if (params.fd <= 0.)
        throw std::runtime_error("Error in BaseGenerator::configure function." + 
                                 std::string(" Invalid parameters fd: ") +
                                 std::to_string(params.fd));

    if (params.f <= 0.)
        throw std::runtime_error("Error in BaseGenerator::configure function." + 
                                 std::string(" Invalid parameters f: ") +
                                 std::to_string(params.f));

    if (params.n == 0)
        throw std::runtime_error("Error in BaseGenerator::configure function." + 
                                 std::string(" Invalid parameters numBits: ") +
                                 std::to_string(params.n));

    if (params.vel <= 0.)
        throw std::runtime_error("Error in BaseGenerator::configure function." + 
                                 std::string(" Invalid parameters infoVel: ") +
                                 std::to_string(params.vel));

    // Koef for normal trasmission (F_info / f << 1)
    double koeff = 1. / params.vel / params.f;
    std::cout << "config. koeff: " << koeff << std::endl;

    if (koeff >= 0.1)
        throw std::runtime_error("Error in BaseGenerator::configure function." + 
                                 std::string(" Invalid parameters infoVel. F_info / f_carrier =  ") +
                                 std::to_string(koeff) + std::string(", (koeff << 1)"));
    if (params.f * 2 > params.fd)
        throw std::runtime_error("Error in BaseGenerator::configure function." + 
                                 std::string(" Invalid parameters fd and f. fd >= 2 * f"));
 
    m_NumBits = params.n;
    m_Type    = params.type;
    
    // T = 1 / fd
    // t = numBits * infoVel
    // Num samples = t / T
    m_NumSampl = params.n * params.vel * params.fd;
    std::cout << "config. numSamples: " << m_NumSampl << std::endl;

    // Num samples / numBits
    m_SamplPerBit = params.vel  * params.fd;
    std::cout << "config. Samples per bit: " << m_SamplPerBit << std::endl;

    // T = 1 / fd
    // phase = ph0 + f * t, where t = n * T
    m_DPhase = params.f / params.fd;
    std::cout << "config. d Phase: " << m_DPhase << std::endl;

    // d_f = Fd / 2, whete Fd - Info freq
    double fMin = 1. / 2. / params.vel;
    std::cout << "d_f: " << fMin << std::endl;

    // m_DPhaseFreqMod = {(params.f + fMin) / params.fd, (params.f - fMin) / params.fd};
    m_DPhaseFreqMod = {(params.f * (1 + 0.5)) / params.fd, (params.f * (1 - 0.5)) / params.fd};

    return;
};

// Generate data signal
//! [out] bits     - Generated input bits
//! [out] data_out - Generated sample data
void BaseGenerator::generate(std::vector<uint8_t>& bits, std::vector<std::complex<double>>& data_out);
{
    // Generate random bits
    m_Gen.generateBits(bits, m_NumBits);

    // Code bits with Gold Seq
    auto codeBits = genCodeBits(bits);

    data_out.resize(codeBits.size() * m_numSamplesPerBit);

    // Generate IQ
    generateIQ(codeBits, data_out);

    return;
}

// Generate Inpulse Responce
void BaseGenerator::generateInpulseResp(std::vector<std::vector<std::complex<double>>>& impulse_res)
{
    impulse_res.resize(4, std::vector<std::complex<double>>(m_SamplPerBit * m_goldSeq.symb1.size()));

    generateIQ(m_goldSeq.symb1, impulse_res[0]);
    generateIQ(m_goldSeq.symb2, impulse_res[1]);
    generateIQ(m_goldSeq.symb3, impulse_res[2]);
    generateIQ(m_goldSeq.symb4, impulse_res[3]);

    return;
}

// Generate IQ
void BaseGenerator::generateIQ(const std::vector<uint8_t>& bits, std::vector<std::complex<double>>& data_out)
{
    // Pointer to current info bit val
    auto cur_bit = bits.begin();

    double I = 0;
    double Q = 0;

    static const uint32_t incrementBits = 2;

    uint32_t size = data_out.size();

    for (uint32_t i = 0; i < size; ++i)
    {
        I = (*cur_bit ? 1 : -1);
        Q = (*(cur_bit + 1) ? 1 : -1);

        data_out[i] = I + Q;

        // We need to change bit val
        if (!(i % m_SamplPerBit) && (i != 0))
            cur_bit += incrementBits;
    }

    return;
}

// Here we assumed what 1 symbol is 2 bits
// We need to get this and replace with Gold Seq  

// Generate Code bits with Gold Seq
std::vector<uint8_t> BaseGenerator::generateCodeBits(const std::vector<uint8_t>& bits)
{
    uint32_t size     = bits.size();
    uint32_t goldSize = m_goldSeq.symb1.size();

    std::vector<uint8_t> res(size * goldSize);
    auto resIterator = res.begin();

    for (uint32_t i = 0; i < size; i += 2)
    {
        if (bits[i] == 0 && bits[i + 1] == 0)          // 00
            res.insert(resIterator, m_goldSeq.symb1.begin(), m_goldSeq.symb1.end());
        else if (bits[i] == 0 && bits[i + 1] == 1)     // 01
            res.insert(resIterator, m_goldSeq.symb2.begin(), m_goldSeq.symb2.end());
        else if (bits[i] == 1 && bits[i + 1] == 0)     // 10
            res.insert(resIterator, m_goldSeq.symb3.begin(), m_goldSeq.symb3.end());
        else                                           // 11
            res.insert(resIterator, m_goldSeq.symb4.begin(), m_goldSeq.symb4.end());

        resIterator += goldSize;
    }

    return res;
}


// Get number of samples per bit
uint32_t BaseGenerator::getNumSamplesPerBit()
{
    return m_SamplPerBit;
}

// Configure Data Processor
void DataProcessor::config(const cfg& params)
{
    m_Cfg = params;
    m_GenData.configure(params);

    if (params.type == SignalType::amplitude)
        fileName = "../data/ber_am.txt";
    else if (params.type == SignalType::phase)
        fileName = "../data/ber_pm.txt";
    else if (params.type == SignalType::freq)
        fileName = "../data/ber_fm.txt";
    else
        throw std::runtime_error("Error in DataProcessor::config! Signal type is ndf!");

    return;
}

// Run Data Processing
void DataProcessor::run(uint32_t num_runs)
{
    size_t counter = 0;
    Correlator corr;

    // Temp data for processing
    std::vector<std::complex<double>> firstSignal;
    std::vector<std::complex<double>> secondSignal;
    std::vector<double>               correlation;
    double   shifted_size_per    = 0.3; // 30 %
    uint32_t max_metric_id       = 0;
    uint32_t shifted_signal_size = 0;
    double persent = 0;

    double dt = 0;

    // Processing steps
    for (uint32_t i = 0; i < num_runs; ++i)
    {
        // Generate large part
        m_GenData.generate(firstSignal);

        shifted_signal_size = shifted_size_per * firstSignal.size();

        // Generate min part
        dt = SignalGenerator::generateShiftedSignal(m_Cfg.fd, shifted_signal_size, firstSignal, m_UniGen.generate(), secondSignal);

        m_Noise.addNoise(firstSignal,  m_Cfg.snr1);
        m_Noise.addNoise(secondSignal, m_Cfg.snr2);

        corr.correlate(firstSignal, secondSignal, correlation, max_metric_id);

        if (std::abs(max_metric_id - dt) < m_GenData.getNumSamplesPerBit())
            counter++;
    }

    persent = (double)counter / (double)num_runs;

    Utils::writeBer(persent, fileName);

    return;
}

// Run Data Processing with writing temp data 
void DataProcessor::run()
{
    std::vector<std::complex<double>> signalIQ;
    std::vector<uint8_t>              bits;
    std::vector<std::vector<double>>  correlation(4);
    std::vector<std::vector<double>>  impulceResponce(4);

    // Generate signal
    m_GenData.generate(bits, signalIQ);

    // Add noise
    m_Noise.addNoise(signalIQ, m_Cfg.snr);

    // Generate Impulse Responce

    Correlator corr;
    auto outBits = corr.correlate(goldSeq, singalIQ, correlation);
    
    Utils::write(bits, std::string("../data/in_bits.txt"));
    Utils::write(outBits, std::string("../data/out_bits.txt"));
    Utils::write(signalIQ, std::string("../data/signalIQ.txt"));

    Utils::write(correlation[0], std::string("../data/correlation1.txt"));
    Utils::write(correlation[1], std::string("../data/correlation2.txt"));
    Utils::write(correlation[2], std::string("../data/correlation3.txt"));
    Utils::write(correlation[3], std::string("../data/correlation4.txt"));

    return;
}

//! M-seq generator functions

//! Get seq
std::vector<uint8_t> MSeqGenerator::getSequence()
{
    // Size: max Size: 2^n - 1. We need (size % 2 == 0)
    std::vector<uint8_t> res((1U << m_N) - 2);

    for (auto& item : res)
        item = getNext();

    return res;
}

//! Get next bit
uint8_t MSeqGenerator::getNext()
{
    uint8_t res = m_Seed & 1U;

    uint32_t feedBack = 0;
    uint32_t regResponce = m_Seed & m_Poly;

    // Get XOR of all responce
    while (regResponce)
    {
        feedBack ^= (regResponce & 1U);
        regResponce >>= 1;
    }

    // Reg shift
    m_Seed >>= 1;

    // New val
    if (feedBack)
        m_Seed |= (1U << (m_N - 1));

    return res;
}


//! Gold Generator functions
// Get seq
GoldSeq GoldGenerator::getSequence()
{
    GoldSeq res;

    // Gen 2 Gold Seq
    auto m1 = m_Gen1->getSequence();
    auto m2 = m_Gen2->getSequence();

    res.symb1 = xorVec(m1, m2);
    shiftVec(m2);
    res.symb2 = xorVec(m1, m2);
    shiftVec(m2);
    res.symb3 = xorVec(m1, m2);
    shiftVec(m2);
    res.symb4 = xorVec(m1, m2);

    return res;
}

//! Shift Vector function
void GoldGenerator::shiftVec(std::vector<uint8_t>& vec)
{
    auto frontEl = vec.front();

    vec.erase(vec.begin());

    vec.push_back(frontEl);

    return;
}

//! Xor vector fucntion
std::vector<uint8_t> GoldGenerator::xorVec(const std::vector<uint8_t>& vec1, const std::vector<uint8_t>& vec2)
{
    auto res = vec1;

    for (uint32_t i = 0; i < res.size(); ++i)
        res[i] ^= vec2[i];

    return res;
}

