cmake_minimum_required(VERSION 3.15)
project(data_processing)

# Установка типа сборки
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# Флаги для отладки
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

# Поиск FFTW
find_path(FFTW3_INCLUDE_DIR fftw3.h)
find_library(FFTW3_LIBRARY NAMES fftw3 fftw3-3)

if (FFTW3_INCLUDE_DIR AND FFTW3_LIBRARY)
    message(STATUS "Found FFTW3: ${FFTW3_LIBRARY}")
else()
    message(FATAL_ERROR "FFTW3 not found!")
endif()

add_executable(${PROJECT_NAME} 
    main.cpp
    Generator.cpp
    Correlator.cpp
)

target_compile_options(${PROJECT_NAME} PRIVATE -pthread)
target_link_libraries(${PROJECT_NAME} PRIVATE ${FFTW3_LIBRARY} pthread)

target_include_directories(${PROJECT_NAME} PRIVATE ${FFTW3_INCLUDE_DIR})

# Убедитесь, что символы отладки включены
target_compile_options(${PROJECT_NAME} PRIVATE $<$<CONFIG:Debug>:-g>)