#pragma once
#include <iostream>
#include <time.h>

#include "Matrix.h"

class Timer {
private:
	clock_t startTime;
	clock_t endTime;
public:
	Timer();
	void start();
	void stop();
	double getTimeMilliSec();
};

class TestMatrix {
private:
	bool compare(Matrix<double>, Matrix<double>);
	void rnd(Matrix<double>*);
public:
	bool testMath(long long);
	double* determinantTime(int, long long);
	bool determinantAccuracy(int, long long);
};

