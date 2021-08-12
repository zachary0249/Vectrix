#ifndef vx_H
#define vx_H

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <cstdlib>
#include <map>

template <typename T>
class vx
{
    // defining class members
private:
    std::vector<T> vectrix;
    unsigned rows;
    unsigned cols;

public:
    vx(unsigned _rows, unsigned _cols, const T &_initial); // takes # rows, # cols, and the inital value for each cell
    vx(const vx<T> &elem);                                 // copy constructor
    vx(const std::vector<T> &elem);
    vx(const std::vector<std::vector<T>> &elem); // takes 2D vector and returns as matrix
    virtual ~vx();                               // deconstructor

    // operator overloading for math operations

    // matrix math
    vx<T> &operator=(const vx<T> &elem);
    vx<T> operator+(const vx<T> &elem);
    vx<T> &operator+=(const vx<T> &elem);
    vx<T> operator*(const vx<T> &elem);
    vx<T> &operator*=(const vx<T> &elem);
    vx<T> operator-(const vx<T> &elem);
    vx<T> &operator-=(const vx<T> &elem);
    vx<T> matmul(const vx<T> &elem1, const vx<T> &elem2);
    vx<T> transpose();
    float sum();
    float mean();
    float median() const;
    float SD();
    std::vector<T> flatten(); // flattens the matrix into a vector
    T max();
    T min();
    vx<T> &absolute();                                                                 // absolute value of a matrix
    std::vector<vx<T>> balance(const vx<T> &elem1, const vx<T> &elem2, const T value); // balances out two matrices such that their result is of the same order
    std::vector<T> diag_vec();                                                         // returns the diagonal of a matrix

    // vector operations
    float dot_product(const vx<T> &elem);   // dot product of two vectors
    float cross_product(const vx<T> &elem); // cross product of two vectors

    // scalar to matrix math
    vx<T> operator+(const T elem);
    vx<T> &operator+=(const T elem);
    vx<T> operator-(const T elem);
    vx<T> &operator-=(const T elem);
    vx<T> operator*(const T elem);
    vx<T> &operator*=(const T elem);
    vx<T> operator/(const T elem);
    vx<T> &operator/=(const T elem);
    vx<T> operator^(const T elem);

    // utility functions
    std::vector<unsigned> dims(bool print = true) const; // for matrix
    unsigned nitems() const;                             // for matrix
    unsigned nrows() const;                              // for matrix
    unsigned ncols() const;                              // for matrix
    vx<T> reshape(const unsigned rows, const unsigned cols);
    bool is_ragged();                                                        // ragged meaning rows of different lengths
    void print() const;                                                      // prints the matrix
    T &operator()(const unsigned &row, const unsigned &col = 0);             // accessing specific index of matrix
    const T &operator()(const unsigned &row, const unsigned &col = 0) const; // accessing specific index of (const) matrix

    unsigned digit_count() const; // max digits of elements in matrix
};

/*
DEFINING NON CLASS FUNCTIONS AND VARIABLES:
*/

static const double pi = 3.14159265358979323846;
static const double e = 2.71828182845904523536;

// gets the index given the cols and rows
unsigned GetIndex(const unsigned &num_cols, const unsigned &_row, const unsigned &_col);

// random number generation
std::vector<double> getRandom_decimal(const unsigned n);

// vector math
template <typename T>
std::vector<T> &add(std::vector<T> &elem1, std::vector<T> &elem2); // add two vectors

template <typename T>
std::vector<T> &add(std::vector<T> &elem1, const T &elem2); // add each member of a vector with a scalar

template <typename T>
std::vector<T> &subtract(std::vector<T> &elem1, const T &elem2); // subtract each member of a vector with a scalar

template <typename T>
std::vector<T> &subtract(std::vector<T> &elem1, std::vector<T> &elem2); // subtract two vectors

template <typename T>
std::vector<T> &scale(std::vector<T> &elem1, const T &elem2); // multiplies each member by the scale factor; to shrink the vector apply the inverse of the scalar

template <typename T>
float sum(const std::vector<T> &elem);

template <typename T>
float sum(const std::vector<std::vector<T>> &elem);

template <typename T>
float mean(const std::vector<T> &elem);

template <typename T>
float mean(const std::vector<std::vector<T>> &elem);

template <typename T>
float SD(std::vector<T> &elem);

template <typename T>
float SD(std::vector<std::vector<T>> &elem);

template <typename T>
vx<T> expand(const std::vector<T> &elem); // expands a 1D vector into a vx row Matrix

template <typename T>
T max(const std::vector<T> &elem);

template <typename T>
T max(const std::vector<std::vector<T>> &elem);

template <typename T>
T min(const std::vector<T> &elem);

template <typename T>
T min(const std::vector<std::vector<T>> &elem);

template <typename T>
float angle_between(const std::vector<T> &elem1, const std::vector<T> &elem2); // angle between two vectors; in degrees

template <typename T>
float distance(const std::vector<T> &elem); // Calculates the distance between two points (endpoints of the line); AKA vector magnitude

template <typename T>
float vec_angle(const std::vector<T> &elem); // gets the angle (direction) of a vector; in degrees

template <typename T>
float slope(const std::vector<T> &elem);

template <typename T>
float rad_to_deg(T value);

template <typename T>
float deg_to_rad(T value);

template <typename T>
std::vector<unsigned> dims(const std::vector<std::vector<T>> &elem, bool print = true); // for 2D vector

template <typename T>
bool is_ragged(const std::vector<std::vector<T>> &elem); // ragged meaning rows of different lengths

template <typename T>
void print(const T &elem); // prints scalar

template <typename T>
void print(const std::vector<T> &elem); // prints vector

template <typename T>
void print(const std::vector<std::vector<T>> &elem); // prints 2D vector

template <typename T>
T &absolute(T &value);

template <typename T>
std::vector<T> &absolute(std::vector<T> &elem);

template <typename T>
std::vector<std::vector<T>> &absolute(std::vector<std::vector<T>> &elem);

template <typename T>
unsigned nitems(const std::vector<std::vector<T>> &elem);

template <typename T>
unsigned nrows(const std::vector<std::vector<T>> &elem);

template <typename T>
unsigned ncols(const std::vector<std::vector<T>> &elem);

template <typename T>
std::vector<T> range1D(const T stop, const T start, const T increment); // range vector

template <typename T>
vx<T> range2D(const T stop, const T start, const T increment, const unsigned numrows); // range matrix

template <typename T>
unsigned digit_count(const T &num);

template <typename T>
unsigned digit_count(const std::vector<T> &elem); // returns the max digit count in vector

template <typename T>
std::vector<T> &sort(std::vector<T> &elem);

template <typename T>
float median(std::vector<T> &elem);

#include "vectrix.cpp"

#endif