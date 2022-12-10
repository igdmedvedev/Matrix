#include <iostream>

#include "Utils.h"
#include "Matrix.h"

int main() {
	srand(time(NULL));
	setlocale(LC_ALL, "Russian");

	TestMatrix test;
	
	std::cout << test.determinantAccuracy(9, 100) << std::endl;

	double* detTime = test.determinantTime(9, 10);
	std::cout << detTime[0] << ' ' << detTime[1] << ' ' << std::endl;

	try {
		test.testMath(1000);
	}
	catch (MyException& ex) {
		std::cout << ex.what();
	}
	
	Matrix<double> m1(3, 3);
	m1(0, 0) = 2; m1(0, 1) = 5; m1(0, 2) = 7;
	m1(1, 0) = 6; m1(1, 1) = 3; m1(1, 2) = 4;
	m1(2, 0) = 5; m1(2, 1) = -2; m1(2, 2) = -3;
	std::cout << m1.inverse();
}



