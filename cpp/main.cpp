#include <iostream>
#include <thread>
#include <vector>
#include <memory>

#include "Generator.h"
#include "Correlator.h"

// Parse Cfg for the demonstration
void parceCfg(cfg& cfg, char* argv[])
{
    cfg.fd   = std::stod(argv[1]);
    cfg.f    = std::stod(argv[2]);
    cfg.n    = std::stoi(argv[3]);
    cfg.vel  = std::stod(argv[4]);
    cfg.dt   = std::stod(argv[5]);
    cfg.snr1 = std::stod(argv[6]);
    cfg.snr2 = std::stod(argv[7]);
    cfg.type = static_cast<SignalType>(std::stoi(argv[8]));
}

// Parse Cfg for the researching
void parceCfgResearch(cfg& cfg, uint32_t& numRans, char* argv[])
{
    cfg.fd   = std::stod(argv[1]);
    cfg.f    = std::stod(argv[2]);
    cfg.n    = std::stoi(argv[3]);
    cfg.vel  = std::stod(argv[4]);
    cfg.snr1 = std::stod(argv[5]);
    cfg.snr2 = std::stod(argv[6]);
    numRans  = std::stoi(argv[7]);
    cfg.is_random_dt = true;
}

// Функция для запуска обработки в потоке
void runProcessor(DataProcessor& processor, const cfg& config, uint32_t numRans)
{
    processor.config(config);
    processor.run(numRans);
}

int main(int argc, char* argv[])
{
    if (argc != 9 && argc != 8)
    {
        std::cerr << "Incorrect input number of parameters: " << argc << std::endl;
        std::cerr << "Usage for demo: " << argv[0] << " fd f n vel dt snr1 snr2 type" << std::endl;
        std::cerr << "Usage for research: " << argv[0] << " fd f n vel snr1 snr2" << std::endl;
        return 1;
    }

    cfg config;
    uint32_t num_runs = 0;
    
    if (argc == 9)
    {
        parceCfg(config, argv);
    }
    else
    {
        parceCfgResearch(config, num_runs, argv);
    }

    switch (argc)
    {
    case 9:  // Demonstration mode
    {
        DataProcessor proc;
        proc.config(config);
        proc.run();
        break;
    }
    case 8:  // Researching mode
    {
        // Создаем процессоры для каждого типа модуляции
        DataProcessor am;  // AM
        DataProcessor pm;  // BPSK
        DataProcessor fm;  // MFM

        // Конфигурации для каждого типа
        cfg am_cfg = config;
        cfg pm_cfg = config;
        cfg fm_cfg = config;

        am_cfg.type = SignalType::amplitude;
        pm_cfg.type = SignalType::phase;
        fm_cfg.type = SignalType::freq;

        // Создаем потоки
        std::vector<std::thread> threads;
        
        // Запускаем обработку в отдельных потоках
        threads.emplace_back(runProcessor, std::ref(am), am_cfg, num_runs);
        threads.emplace_back(runProcessor, std::ref(pm), pm_cfg, num_runs);
        threads.emplace_back(runProcessor, std::ref(fm), fm_cfg, num_runs);

        // Ожидаем завершения всех потоков
        for (auto& thread : threads)
            if (thread.joinable())
                thread.join();

        std::cout << "All modulation types processed successfully!" << std::endl;
        break;
    }
    default:
        std::cerr << "Unexpected number of arguments" << std::endl;
        return 1;
    }

    return 0;
}