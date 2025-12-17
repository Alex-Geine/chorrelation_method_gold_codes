#include "Generator.h"
#include "Correlator.h"
#include <fstream>
#include <string>
#include <type_traits>
#include <algorithm>

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

    if constexpr (std::is_same_v<T, std::complex<double>>)
    {
        for (uint32_t i = 0; i < data.size() - 1; ++i)
            outFile << data[i].real() << ", " << data[i].imag() << ", ";
        outFile << data[data.size() - 1].real() << ", " << data[data.size() - 1].imag() << " " << std::endl;
    }
    else if constexpr (std::is_same_v<T, uint8_t>)
    {
        for (uint32_t i = 0; i < data.size() - 1; ++i)
            outFile << (uint32_t)data[i] << ", ";
        outFile << (uint32_t)data[data.size() - 1];
    }
    else
    {
        for (uint32_t i = 0; i < data.size() - 1; ++i)
            outFile << data[i] << ", ";
        outFile << data[data.size() - 1] << std::endl;
    }

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

// Add noise in data
//! [in/out] data     - Input/Output data
//! [in]     snr      - Signal to Noise Ratio
void NoiseInjector::addNoise(std::vector<std::complex<double>>& data, double snr)
{
    double energySignal = 0;
    double energyNoise  = 0;

    for (auto& it: data)
        energySignal += std::norm(it);

    // snr = log(energy / noise)
    // 10^snr = enegry / noise
    // noise = energy / 10^snr
    double lin_snr = std::pow(10, snr / 10.);
    double noise_power = 1. / std::sqrt(2.);

    gen.config(noise_power);

    uint32_t size = data.size();

    std::vector<std::complex<double>> noise;
    gen.generateAwgn(noise, size);

    for (auto& it: noise)
        energyNoise += std::norm(it);

    double alfa = std::sqrt(energySignal / energyNoise / lin_snr);
    
    for (uint32_t i = 0; i < size; ++i)
        data[i] += alfa * noise[i];

    energyNoise = 0;
    for (auto& it: noise)
        energyNoise += std::norm(it * alfa);

    std::cout << "SNR[dB]: " << snr << "real SNR[dB]: " << 10. * log10(energySignal / energyNoise) << std::endl;

    // Utils::writeBer(10. * log10(energySignal / energyNoise), "real_snr.txt");

    return;
}

//! Configurate signal generator
//! [in] params - Configuration parameters
//! [in] seq    - Gold Sequence
void BaseGenerator::configure(const cfg& params, const GoldSeq& seq)
{
    if (params.fd <= 0.)
        throw std::runtime_error("Error in BaseGenerator::configure function." + 
                                 std::string(" Invalid parameters fd: ") +
                                 std::to_string(params.fd));

    if (params.n == 0)
        throw std::runtime_error("Error in BaseGenerator::configure function." + 
                                 std::string(" Invalid parameters numBits: ") +
                                 std::to_string(params.n));

    if (params.vel <= 0.)
        throw std::runtime_error("Error in BaseGenerator::configure function." + 
                                 std::string(" Invalid parameters infoVel: ") +
                                 std::to_string(params.vel));

    m_NumBits = params.n;

    // Num samples / numBits
    m_SamplPerSymb = 1. / params.vel  *  params.fd;

    // Get Gold size
    m_goldSeq = seq;

    return;
}

// Generate data signal
//! [out] bits     - Generated input bits
//! [out] data_out - Generated sample data
void BaseGenerator::generate(std::vector<uint8_t>& bits, std::vector<std::complex<double>>& data_out)
{
    // Generate random bits
    m_Gen.generateBits(bits, m_NumBits);

    // Code bits with Gold Seq
    std::vector<uint8_t> codeBits = generateCodeBits(bits);

    data_out.resize(codeBits.size() * m_SamplPerSymb / 2);

    // Generate IQ
    generateIQ(codeBits, data_out);

    std::cout << "codeBits size: " << codeBits.size() << std::endl;
    std::cout << "dataOut size: " << data_out.size() << std::endl;

    return;
}

// Generate Inpulse Responce
void BaseGenerator::generateInpulseResp(std::vector<std::vector<std::complex<double>>>& impulse_res)
{
    impulse_res.resize(4, std::vector<std::complex<double>>(m_SamplPerSymb * m_goldSeq.symb1.size() / 2));

    std::cout << "impulse responce: " << m_SamplPerSymb * m_goldSeq.symb1.size() / 2 << std::endl;

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

        data_out[i] = {I, Q};

        // We need to change bit val
        if (!(i % m_SamplPerSymb) && (i != 0))
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

    std::vector<uint8_t> res(size * goldSize / 2);
    uint32_t resIterator = 0;

    for (uint32_t i = 0; i < size; i += 2)
    {
        const std::vector<uint8_t>* symbol = nullptr;

        if (bits[i] == 0 && bits[i + 1] == 0)
        {
            symbol = &m_goldSeq.symb1;
        }
        else if (bits[i] == 0 && bits[i + 1] == 1)
        {
            symbol = &m_goldSeq.symb2;
        }
        else if (bits[i] == 1 && bits[i + 1] == 0)
        {
            symbol = &m_goldSeq.symb3;
        }
        else {
            symbol = &m_goldSeq.symb4;
        }

        std::copy(symbol->begin(), symbol->end(), res.begin() + resIterator);
        resIterator += goldSize;
    }

    std::cout << "res.size(): " << res.size() << std::endl;

    return res;
}

// Configure Data Processor
void DataProcessor::config(const cfg& params)
{
    m_Cfg = params;


    Utils::writeBer(params.poly1, "poly1.txt");
    Utils::writeBer(params.poly2, "poly2.txt");
    

    std::mt19937 generator(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<uint32_t> distribution(1, 0xFFFFFFFF);
    static const uint32_t s_RegSize = 5;
    GoldGenerator gen(s_RegSize, distribution(generator), distribution(generator), params.poly1, params.poly2);     

    auto goldCode = gen.getSequence();

    m_GenData.configure(params, goldCode);

    fileName = "../data/ber_qpsk.txt";

    return;
}

// Run Data Processing
void DataProcessor::run(uint32_t num_runs)
{
    size_t counter = 0;
    double persent = 0;
    Correlator corr;

    // Temp data for processing
    std::vector<std::complex<double>> signalIQ;
    std::vector<uint8_t>              bits;
    std::vector<uint8_t>              outBits;
    std::vector<std::vector<double>>  correlation(4);
    std::vector<std::vector<std::complex<double>>> impulceResponce;


    Utils::writeBer(num_runs, "num_runs.txt");

    // Generate impulse responce
    m_GenData.generateInpulseResp(impulceResponce);

    // Processing steps
    for (uint32_t i = 0; i < num_runs; ++i)
    {
        // Generate signal
        m_GenData.generate(bits, signalIQ);

        // Add noise
        m_Noise.addNoise(signalIQ, m_Cfg.snr);

        // Correlate
        outBits = corr.correlate(impulceResponce, signalIQ, correlation);

        if (std::equal(bits.begin(), bits.end(), outBits.begin()))
            counter++;
    }

    Utils::writeBer(counter, "counter.txt");


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
    std::vector<std::vector<std::complex<double>>> impulceResponce;

    // Generate signal
    m_GenData.generate(bits, signalIQ);

    std::cout << "bits: " << std::endl;
    for (auto it : bits)
        std::cout << (uint32_t)it << ", ";
    std::cout << std::endl;

    // Generate impulse responce
    m_GenData.generateInpulseResp(impulceResponce);

    // Add noise
    m_Noise.addNoise(signalIQ, m_Cfg.snr);

    // Correlate data
    Correlator corr;
    auto outBits = corr.correlate(impulceResponce, signalIQ, correlation);
    
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
    // Size: max Size: 2^n - 1
    std::vector<uint8_t> res((1U << m_N) - 1);

    for (auto& item : res)
        item = getNext();

    return res;
}

//! Get next bit
uint8_t MSeqGenerator::getNext()
{
   uint8_t res = m_Seed & 1U;
    
    // More efficient feedback calculation using popcount
    uint32_t feedbackBits = m_Seed & m_Poly;
    
    // Count set bits (parity)
    uint32_t feedback = 0;
    while (feedbackBits) {
        feedback ^= 1;
        feedbackBits &= feedbackBits - 1;  // Clear lowest set bit
    }
    
    // Shift register
    m_Seed >>= 1;
    
    // Insert feedback at MSB position
    if (feedback) {
        m_Seed |= (1U << (m_N - 1));
    }

    return res;
}


//! Gold Generator functions
// Get seq
GoldSeq GoldGenerator::getSequence()
{
    GoldSeq res;

    // Generate base M-sequences
    auto m1 = m_Gen1->getSequence();
    auto m2 = m_Gen2->getSequence();

    res.symb1 = xorVec(m1, m2);
    
    auto m2_shift1 = m2;
    shiftVec(m2_shift1, 1);
    res.symb2 = xorVec(m1, m2_shift1);
    
    auto m2_shift2 = m2;
    shiftVec(m2_shift2, 2);
    res.symb3 = xorVec(m1, m2_shift2);
    
    auto m2_shift3 = m2;
    shiftVec(m2_shift3, 3);
    res.symb4 = xorVec(m1, m2_shift3);

    res.symb1.erase(res.symb1.end() - 1);
    res.symb2.erase(res.symb2.end() - 1);
    res.symb3.erase(res.symb3.end() - 1);
    res.symb4.erase(res.symb4.end() - 1);

    std::cout << "seq 1: " << std::endl;
    for (auto it : res.symb1)
        std::cout << (uint32_t)it << ", ";
    std::cout << std::endl;

    std::cout << "seq 2: " << std::endl;
    for (auto it : res.symb2)
        std::cout << (uint32_t)it << ", ";
    std::cout << std::endl;

    std::cout << "seq 3: " << std::endl;
    for (auto it : res.symb3)
        std::cout << (uint32_t)it << ", ";
    std::cout << std::endl;

    std::cout << "seq 4: " << std::endl;
    for (auto it : res.symb4)
        std::cout << (uint32_t)it << ", ";
    std::cout << std::endl;


    return res;
}

//! Shift Vector function
void GoldGenerator::shiftVec(std::vector<uint8_t>& vec, uint32_t shift)
{
    shift %= vec.size();

    std::rotate(vec.begin(), vec.begin() + shift, vec.end());

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

