#pragma once
#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <iomanip>
#include <math.h>
class MyException : public std::exception {
public:
	MyException(const char* msg) : exception(msg) {}
};

template <class T> class Matrix {
private:
	T** matrix;
	size_t rows, cols;

	void clear();
	T*& operator[](size_t);
public:
	Matrix(size_t, size_t, const T& = T());
	Matrix(const Matrix&);
	~Matrix();

	size_t getRows() { return rows; }
	size_t getCols() { return cols; }

	friend std::ostream& operator << (std::ostream& os, const Matrix<T>& obj) {
		for (size_t i = 0; i < obj.rows; i++) {
			for (size_t j = 0; j < obj.cols; j++) {
				std::cout << std::setw(3) << obj.matrix[i][j] << ' ';
			}
			std::cout << std::endl;
		}
		return os;
	}

	T& operator()(size_t, size_t);
	
	Matrix operator*(const T);
	Matrix operator/(const T);
	Matrix operator*(const Matrix&);
	Matrix operator+(const Matrix&);
	Matrix operator-(const Matrix&);
	Matrix operator*=(const T);
	Matrix operator/=(const T);
	Matrix operator*=(const Matrix&);
	Matrix operator+=(const Matrix&);
	Matrix operator-=(const Matrix&);
	Matrix operator=(const Matrix&);

	Matrix transposition();

	double determinant();
	double determinantAlgebraic();
	double determinantGauss();

	double minor(size_t, size_t);
	double algebraicComplement(size_t, size_t);
	Matrix<double> inverse();
};

template <typename T> Matrix<T>::Matrix(size_t rows, size_t cols, const T& val) {
	if (cols < 1 || rows < 1) {
		MyException exception("Ñonstructor failed, invalid indices");
		throw exception;
	}
	this->rows = rows;
	this->cols = cols;

	matrix = new T * [rows];
	for (size_t i = 0; i < rows; i++) {
		matrix[i] = new T[cols];

		for (size_t j = 0; j < cols; j++) {
			matrix[i][j] = val;
		}
	}
}
template <typename T> Matrix<T>::Matrix(const Matrix& obj) {
	rows = obj.rows;
	cols = obj.cols;

	matrix = new T * [rows];
	for (size_t i = 0; i < rows; i++) {
		matrix[i] = new T[cols];

		for (size_t j = 0; j < cols; j++) {
			matrix[i][j] = obj.matrix[i][j];
		}
	}
}
template <typename T> Matrix<T>::~Matrix() {
	clear();
}

template <typename T> void Matrix<T>::clear() {
	for (size_t i = 0; i < rows; i++) {
		delete[] matrix[i];
	}
	delete[] matrix;
}

template <typename T> T& Matrix<T>::operator()(size_t rows, size_t cols) {
	if (rows < 0 || rows > this->rows) {
		MyException exception("Invalid line number");
		throw exception;
	}
	if (cols < 0 || cols > this->cols) {
		MyException exception("Invalid column number");
		throw exception;
	}
	return matrix[rows][cols];
}
template <typename T> T*& Matrix<T>::operator[](size_t n) {
	return matrix[n];
}
template <typename T> Matrix<T> Matrix<T>::operator*(const T num) {
	Matrix<T> tmp = *this;
	tmp *= num;
	return tmp;
}
template <typename T> Matrix<T> Matrix<T>::operator/(const T num) {
	Matrix<T> tmp = *this;
	tmp /= num;
	return tmp;
}
template <typename T> Matrix<T> Matrix<T>::operator*(const Matrix& obj) {
	Matrix<T> tmp = *this;
	tmp *= obj;
	return tmp;
}
template <typename T> Matrix<T> Matrix<T>::operator+(const Matrix& obj) {
	Matrix<T> tmp = *this;
	tmp += obj;
	return tmp;
}
template <typename T> Matrix<T> Matrix<T>::operator-(const Matrix& obj) {
	Matrix<T> tmp = *this;
	tmp -= obj;
	return tmp;
}
template <typename T> Matrix<T> Matrix<T>::operator*=(const T num) {
	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < cols; j++) {
			matrix[i][j] *= num;
		}
	}
	return *this;
}
template <typename T> Matrix<T> Matrix<T>::operator/=(const T num) {
	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < cols; j++) {
			matrix[i][j] /= num;
		}
	}
	return *this;
}
template <typename T> Matrix<T> Matrix<T>::operator*=(const Matrix& obj) {
	if (cols != obj.rows) {
		MyException exception("Unsuitable dimensions");
		throw exception;
	}
	Matrix<T> tmp(rows, obj.cols);

	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < cols; j++) {
			for (size_t k = 0; k < obj.rows; k++) {
				tmp.matrix[i][j] += matrix[i][k] * obj.matrix[k][j];
			}
		}
	}
	*this = tmp;
	return *this;
}
template <typename T> Matrix<T> Matrix<T>::operator+=(const Matrix& obj) {
	if (rows != obj.rows) {
		MyException exception("Unequals strings");
		throw exception;
	}
	if (cols != obj.cols) {
		MyException exception("Unequals columns");
		throw exception;
	}
	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < cols; j++) {
			matrix[i][j] += obj.matrix[i][j];
		}
	}
	return *this;
}
template <typename T> Matrix<T> Matrix<T>::operator-=(const Matrix& obj) {
	if (rows != obj.rows) {
		MyException exception("Unequals strings");
		throw exception;
	}
	if (cols != obj.cols) {
		MyException exception("Unequals columns");
		throw exception;
	}
	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < cols; j++) {
			matrix[i][j] -= obj.matrix[i][j];
		}
	}
	return *this;
}
template <typename T> Matrix<T> Matrix<T>::operator=(const Matrix& obj) {
	clear();
	new (this) Matrix(obj);
	return *this;
}

template <typename T> Matrix<T> Matrix<T>::transposition() {
	Matrix<T> tmp(cols, rows);
	for (size_t i = 0; i < cols; i++) {
		for (size_t j = 0; j < rows; j++) {
			tmp.matrix[i][j] = matrix[j][i];
		}
	}
	return tmp;
}

template <typename T> double Matrix<T>::determinant() {
	if (rows == 1 || rows == 2 || rows >= 10) {
		return determinantAlgebraic();
	}
	else {
		return determinantGauss();
	}
}
template <typename T> double Matrix<T>::determinantAlgebraic() {
	if (rows != cols) {
		MyException exception("Non-square matrix");
		throw exception;
	}

	if (rows == 1) {
		return matrix[0][0];
	}
	else if (rows == 2) {
		return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
	}
	else {
		double d = 0;
		for (size_t k = 0; k < rows; k++) {
			Matrix<T> b(rows - 1, rows - 1);
			for (size_t i = 1; i < rows; i++) {
				for (size_t j = 0; j < rows; j++) {
					if (j == k) {
						continue;
					}
					else if (j < k) {
						b(i - 1, j) = matrix[i][j];
					}
					else {
						b(i - 1, j - 1) = matrix[i][j];
					}
				}
			}
			d += (k % 2 == 0 ? 1 : -1) * matrix[0][k] * b.determinantAlgebraic();
		}
		return d;
	}
}
template <typename T> double Matrix<T>::determinantGauss() {
	if (rows != cols) {
		MyException exception("Non-square matrix");
		throw exception;
	}

	const double EPS = 1e-9;
	double det = 1;

	Matrix<double> tmp(rows, cols);
	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < rows; j++) {
			tmp(i, j) = matrix[i][j];
		}
	}

	for (size_t i = 0; i < rows; ++i) {
		size_t k = i;
		for (size_t j = i + 1; j < rows; ++j) {
			if (fabs(tmp(j, i)) > fabs(tmp(k, i))) {
				k = j;
			}
		}
		if (fabs(tmp(k, i)) < EPS) {
			det = 0;
			break;
		}
		std::swap(tmp[i], tmp[k]);
		if (i != k) {
			det = -det;
		}
		det *= tmp(i, i);
		for (size_t j = i + 1; j < rows; ++j) {
			tmp(i, j) /= tmp(i, i);
		}
		for (size_t j = 0; j < rows; ++j) {
			if (j != i && fabs(tmp(j, i)) > EPS) {
				for (size_t k = i + 1; k < rows; ++k) {
					tmp(j, k) -= tmp(i, k) * tmp(j, i);
				}
			}
		}
	}
	return det;
}

template <typename T> double Matrix<T>::minor(size_t rows, size_t cols) {
	if (this->rows != this->cols) {
		MyException exception("Non-square matrix");
		throw exception;
	}

	size_t n = this->rows - 1;
	Matrix<T> tmp(n, n);

	for (size_t i = 0; i < n + 1; i++) {
		if (i == rows) {
			continue;
		}
		for (size_t j = 0; j < n + 1; j++) {
			if (j == cols) {
				continue;
			}
			if (i < rows && j < cols) {
				tmp(i, j) = matrix[i][j];
			}
			else if (i < rows && j > cols) {
				tmp(i, j - 1) = matrix[i][j];
			}
			else if (i > rows && j < cols) {
				tmp(i - 1, j) = matrix[i][j];
			}
			else {
				tmp(i - 1, j - 1) = matrix[i][j];
			}
		}
	}
	return tmp.determinant();
}
template <typename T> double Matrix<T>::algebraicComplement(size_t rows, size_t cols) {
	return ((cols + rows) % 2 == 0 ? 1 : -1) * minor(rows, cols);
}
template <typename T> Matrix<double> Matrix<T>::inverse() {
	if (rows != cols) {
		MyException exception("Non-square matrix");
		throw exception;
	}
	if (determinant() == 0) {
		MyException exception("Degenerate matrix");
		throw exception;
	}

	Matrix<double> tmp(rows, cols);
	if (rows && cols == 1) {
		tmp(0, 0) = 1 / matrix[0][0];
		return tmp;
	}

	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < cols; j++) {
			tmp(i, j) = algebraicComplement(i, j);
		}
	}
	return tmp.transposition() / determinant();
}

#endif