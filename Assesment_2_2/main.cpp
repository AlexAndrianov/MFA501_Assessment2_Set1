#include <iostream>
#include <fstream>

#include "Equation.h"

using namespace std;

int main()
{
    cout << "Program will read the integral limits and tolerance from args.txt" << endl;
    cout << "Please insert angle values in radians!!!" << endl;

    std::ifstream inputFile("args.txt");

    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return 1;
    }

    std::string equationLine;
    std::getline(inputFile, equationLine);

    std::string from;
    std::getline(inputFile, from);

    std::string to;
    std::getline(inputFile, to);

    std::string tol;
    std::getline(inputFile, tol);

    math::Equation equation;
    equation.parse(equationLine);
    cout << "Parsed equation: " << equation._sintaxis_tree_root->toString() << std::endl;

    double fromd = std::stod(from);
    double tod = std::stod(to);
    double told = std::stod(tol);

    cout << "Integral summ from: " << std::to_string(fromd) << " to: " << std::to_string(tod) << " and tolerance: " << std::to_string(tod) << std::endl;
    cout << equation.calculateIntegral(fromd, tod, told) << std::endl;

    return 0;
}
