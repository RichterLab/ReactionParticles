#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>

#include <cxxopts.hpp>

int main( int argc, char* argv[] ) {
    std::string sVelocity;
    unsigned int FieldWidth, FieldHeight;
    std::vector<std::vector<double>> U, V;

    cxxopts::Options options(argv[0]);
    options.add_options()
        ("v,velocity", "Velocity file", cxxopts::value<std::string>(sVelocity)
            ->default_value("vel1.mat")->implicit_value("vel1.mat"))
        ("w,width", "Velocity field width", cxxopts::value<unsigned int>(FieldWidth)
            ->default_value("50")->implicit_value("50"))
        ("h,height", "Velocity field height", cxxopts::value<unsigned int>(FieldHeight)
            ->default_value("5")->implicit_value("5"))
        ("help", "Print help");
    options.parse(argc, argv);

    if (options.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    // Resize Velocity Fields
    U.resize(FieldHeight);
    for( size_t i = 0; i < FieldHeight; i++ ){
        U[i].resize(FieldWidth);
    }

    V.resize(FieldHeight);
    for( size_t i = 0; i < FieldHeight; i++ ){
        V[i].resize(FieldWidth);
    }

    // Parse Velocity File
    std::ifstream fVelocity(sVelocity, std::ifstream::in);
    for( size_t i = 0; i < FieldHeight; i++ ){
        for( size_t j = 0; j < FieldWidth; j++ ){
            fVelocity >> U[i][j];
        }
    }

    for( size_t i = 0; i < FieldHeight; i++ ){
        for( size_t j = 0; j < FieldWidth; j++ ){
            fVelocity >> V[i][j];
        }
    }

    fVelocity.close();

    std::cout << "U" << std::endl;
    for( size_t i = 0; i < FieldHeight; i++ ){
        for( size_t j = 0; j < FieldWidth; j++ ){
            std::cout << U[i][j] << " ";
        }
        std::cout << std::endl;
    }
}