#ifndef _GENERATOR_H_
#define _GENERATOR_H_

#include <iostream>
#include <vector>
#include <random>
#include <complex>
#include <cstdint>
#include <memory>

struct cfg
{
    double   fd    = 0;    // Sample freq
    uint32_t n     = 0;    // Num info bits
    double   vel   = 0;    // Info velocity
    double   snr   = 0;    // SNR for signal 1
    uint32_t poly1 = 0;  // Init Polynom for Gold code
    uint32_t poly2 = 0;  // Init Polynom for Gold code
};

struct GoldSeq
{
    std::vector<uint8_t> symb1;  // 00
    std::vector<uint8_t> symb2;  // 01
    std::vector<uint8_t> symb3;  // 10
    std::vector<uint8_t> symb4;  // 11
};

//! M-seq generator
class MSeqGenerator
{
    private:

    uint32_t m_Seed = 0; //! Initial Reg Seed
    uint32_t m_Poly = 0; //! Polynom
    uint32_t m_N    = 0; //! Register size

    public:
    
    //! Constructor
    MSeqGenerator(uint32_t n, uint32_t seed, uint32_t poly) : m_Seed((seed & ((1U << n) - 1))), m_Poly((poly & ((1U << n) - 1))), m_N(n)
    {
        if ((n == 0) || (seed == 0) || (poly == 0))
        {
            std::cerr << "Error in MSeqGenerator() function. Invalid arguments: n: " << n << ", seed: " << seed << ", poly: " << poly << std::endl;
            abort();
        }
    }

    //! Get seq
    std::vector<uint8_t> getSequence();

    private:

    //! Get next bit
    uint8_t getNext();
};

//! Gold Sequence Generator
class GoldGenerator
{
    private:

    std::unique_ptr<MSeqGenerator> m_Gen1;  // First M-seq gen  
    std::unique_ptr<MSeqGenerator> m_Gen2;  // Second M-seq gen

    public:

    // Constructor
    GoldGenerator(uint32_t n, uint32_t seed1, uint32_t seed2, uint32_t poly1, uint32_t poly2)
    {
        m_Gen1 = std::make_unique<MSeqGenerator>(n, seed1, poly1);
        m_Gen2 = std::make_unique<MSeqGenerator>(n, seed2, poly2);
    }

    // Get seq
    GoldSeq getSequence();

    private:

    //! Shift Vector function
    void shiftVec(std::vector<uint8_t>& vec, uint32_t numShifts = 1);

    //! Xor vector fucntion
    std::vector<uint8_t> xorVec(const std::vector<uint8_t>& vec1, const std::vector<uint8_t>& vec2);
};

class RandomGenerator
{
private: // variables

    std::mt19937                      gen{};    // Random generator
    std::normal_distribution<double>  dist{};   // Normal distribution

public: // functions

    // Default constructor
    RandomGenerator() : gen(std::random_device{}()) {};

    // Configure function
    //! [in] std  - Standart noise deviation
    //! [in] mean - Mean value 
    void config(double std = 1., double mean = 0.);

    // Generate random bits
    //! [out] data_out - Output generated Bits
    //! [in]  size     - Size of the output data
    template <typename Data_t>
    void generateBits(std::vector<Data_t>& data_out, uint32_t size);

    // Generate AWGN
    //! [out] data_out  - Output generated AWGN
    //! [in]  size      - Size of the output data
    void generateAwgn(std::vector<std::complex<double>>& data_out, uint32_t size);

    // Add noise in data
    //! [in/out] data     - Input/Output data
    //! [in]     snr      - Signal to Noise Ratio
    void addNoise(std::vector<std::complex<double>>& data, double snr);
};

// Generate random bits
//! [out] data_out - Output generated Bits
//! [in]  size     - Size of the output data
template <typename Data_t>
void RandomGenerator::generateBits(std::vector<Data_t>& data_out, uint32_t size)
{
    if (size == 0)
        throw std::runtime_error("Error in generateBits function! Unsupported data size: " + size);

    data_out.resize(size);
    for (auto& it:data_out)
        it = ((dist(gen) > 0) ? static_cast<Data_t>(1) : static_cast<Data_t>(0));

    return;
}

class NoiseInjector
{
private:  // variables
    RandomGenerator gen;    //! Random Generator

public:  // functions

    // Add noise in data
    //! [in/out] data     - Input/Output data
    //! [in]     snr      - Signal to Noise Ratio
    void addNoise(std::vector<std::complex<double>>& data, double snr);
};

class BaseGenerator
{
    private: //! variables

    RandomGenerator m_Gen = {};    //! Random generator to create random bits

    uint32_t m_NumBits     = 0;               //! Number of info bits
    uint32_t m_SamplPerSymb = 0;               //! Number samples per bit
    GoldSeq  m_goldSeq     = {};              //! Gold Sequce

    public: //! fucntions

    //! Configurate signal generator
    //! [in] params - Configuration parameters
    //! [in] seq    - Gold Sequence
    void configure(const cfg& params, const GoldSeq& seq);

    // Generate data signal
    //! [out] bits     - Generated input bits
    //! [out] data_out - Generated sample data
    void generate(std::vector<uint8_t>& bits, std::vector<std::complex<double>>& data_out);

    // Generate Inpulse Responce
    void generateInpulseResp(std::vector<std::vector<std::complex<double>>>& impulse_res);

    private:

    // Generate IQ
    void generateIQ(const std::vector<uint8_t>& bits, std::vector<std::complex<double>>& data_out);

    // Generate Code bits with Gold Seq
    std::vector<uint8_t> generateCodeBits(const std::vector<uint8_t>& bits);
};

class DataProcessor
{
    private:

    BaseGenerator                   m_GenData;   // Generator of data
    NoiseInjector                   m_Noise;     // Noise injector
    cfg                             m_Cfg;       // Configuration params
    std::unique_ptr<GoldGenerator>  m_GoldGen;   // Gold seq Generator

    std::string fileName;                        // Filename to write ber data

    public: // functions

    // Configure Data Processor
    void config(const cfg& params);

    // Run Data Processing
    void run(uint32_t num_runs);

    // Run Data Processing with writing temp data 
    void run();
};


#endif //_GENERATOR_H_
