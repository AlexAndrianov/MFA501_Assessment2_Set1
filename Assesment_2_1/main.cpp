#include <iostream>
#include <fstream>
#include <sstream>

#include "Matrix.h"

using namespace std;

math::Matrix readMatrix(std::ifstream &inputFile)
{
    math::Matrix mtx;

    std::string line;
    while (std::getline(inputFile, line)) {
        std::istringstream lineStream(line);

        std::vector<double> values;
        double value;

        while (lineStream >> value) {
            values.push_back(value);
        }

        mtx.addRow(math::Vector(std::move(values)));
    }

    return mtx;
}

int main()
{
    cout << "Program will read the matrix file" << endl;

    std::ifstream inputFile("matrix.txt");

    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return 1;
    }

    auto matrix = readMatrix(inputFile);

    if(matrix.isEmpty())
    {
        std::cerr << "Matrix should not be empty" << std::endl;
        return 1;
    }

    if(!matrix.isSquare())
    {
        std::cerr << "Matrix should be square" << std::endl;
        return 1;
    }

    try {
        auto eigenValues = matrix.eigenValuesByQRAlgorithm();

        std::cout << "Calculated eigen values: " << std::endl;
        std::cout << eigenValues;
    }
    catch(std::exception &ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}
