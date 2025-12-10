#include <iostream>
#include <thread>
#include <vector>
#include <memory>

#include "Generator.h"
#include "Correlator.h"

// Parse Cfg for the demonstration
void parceCfg(cfg& cfg, char* argv[])
{
    cfg.fd     = std::stod(argv[1]);
    cfg.n      = std::stoi(argv[2]);
    cfg.vel    = std::stod(argv[3]);
    cfg.snr    = std::stod(argv[4]);
    cfg.poly1  = std::stoi(argv[5]);
    cfg.poly2  = std::stoi(argv[6]);
}

// Parse Cfg for the researching
void parceCfgResearch(cfg& cfg, uint32_t& numRuns, char* argv[])
{
    cfg.fd     = std::stod(argv[1]);
    cfg.n      = std::stoi(argv[2]);
    cfg.vel    = std::stod(argv[3]);
    cfg.snr    = std::stod(argv[4]);
    cfg.poly1  = std::stoi(argv[5]);
    cfg.poly2  = std::stoi(argv[6]);
    numRuns    = std::stoi(argv[5]);
}

// Функция для запуска обработки в потоке
void runProcessor(DataProcessor& processor, const cfg& config, uint32_t numRuns)
{
    processor.config(config);
    processor.run(numRuns);
}

int main(int argc, char* argv[])
{
    if (argc != 7 && argc != 8)
    {
        std::cerr << "Incorrect input number of parameters: " << argc << std::endl;
        std::cerr << "Usage for demo: " << argv[0] << " fd, n, vel, snr" << std::endl;
        std::cerr << "Usage for research: " << argv[0] << " fd, n, vel, snr, num runs" << std::endl;
        return 1;
    }

    cfg config;
    uint32_t num_runs = 0;
    
    if (argc == 7)
    {
        parceCfg(config, argv);
    }
    else
    {
        parceCfgResearch(config, num_runs, argv);
    }

    switch (argc)
    {
    case 7:  // Demonstration mode
    {
        DataProcessor proc;
        proc.config(config);
        proc.run();
        break;
    }
    case 8:  // Researching mode
    {
        DataProcessor proc;
        proc.config(config);
        proc.run(num_runs);
        break;
    }
    default:
        std::cerr << "Unexpected number of arguments" << std::endl;
        return 1;
    }

    return 0;
}