#include "Utils.h"

Timer::Timer() {
	startTime = -1;
	endTime = -1;
}
void Timer::start() {
	startTime = clock();
	endTime = -1;
}
void Timer::stop() {
	endTime = clock();
}
double Timer::getTimeMilliSec()
{
	return (double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC;
}

bool TestMatrix::compare(Matrix<double> myMatrix1, Matrix<double> myMatrix2) {
	if (myMatrix1.getRows() != myMatrix2.getRows()) {
		return false;
	}
	if (myMatrix1.getCols() != myMatrix2.getCols()) {
		return false;
	}

	for (size_t i = 0; i < myMatrix1.getRows(); i++) {
		for (size_t j = 0; j < myMatrix1.getCols(); j++) {
			if (abs(myMatrix1(i, j) - myMatrix2(i, j)) > 1e-8) {
				return false;
			}
		}
	}

	return true;
}
void TestMatrix::rnd(Matrix<double>* myMatrix) {
	for (size_t i = 0; i < myMatrix->getRows(); i++) {
		for (size_t j = 0; j < myMatrix->getCols(); j++) {
			(*myMatrix)(i, j) = ((double)rand() / RAND_MAX) * 10 - 5;
		}
	}
}

bool TestMatrix::testMath(long long iterations) {
	for (int k = 0; k < iterations; k++) {
		size_t n = rand() % 5 + 1;

		Matrix<double> A(n, n);
		rnd(&A);

		Matrix<double> E(n, n);
		for (size_t i = 0; i < E.getCols(); i++) {
			E(i, i) = 1.0;
		}

		if (!compare(A * A.inverse(), E)) {
			MyException exception("Fail");
			throw exception;
			std::cout << E << std::endl;
		}
	}
	return true;
}
double* TestMatrix::determinantTime(int sizeMatrix, long long iterations) {
	Matrix<double>* m1 = new Matrix<double>(sizeMatrix, sizeMatrix);

	Timer* timer = new Timer();
	double time[2] = { 0 };
	
	for (int i = 0; i < iterations; i++) {
		rnd(m1);

		timer->start();
		m1->determinantAlgebraic();
		timer->stop();
		time[0] += timer->getTimeMilliSec();

		timer->start();
		m1->determinantGauss();
		timer->stop();
		time[1] += timer->getTimeMilliSec();
	}
	
	return time;
}

bool TestMatrix::determinantAccuracy(int sizeMatrix, long long iterations) {
	Matrix<double> m1(sizeMatrix, sizeMatrix);

	for (int i = 0; i < iterations; i++) {
		rnd(&m1);
		double det = m1.determinantAlgebraic();
		double det2 = m1.determinantGauss();
		if (abs(det - det2) > 1e-7) {
			printf("%.10f %.10f\n", det, det2);
			return false;
		}
	}
	return true;
}