#ifndef vx_CPP
#define vx_CPP

#include "vectrix.h"

unsigned GetIndex(const unsigned &num_cols, const unsigned &_row, const unsigned &_col)
{
    if (num_cols == 0 || _col == 0)
    {
        return _row;
    }
    else
    {
        return _row * num_cols + _col;
    }
}

// parameter constructor
template <typename T>
vx<T>::vx(const std::vector<T> &elem)
{
    this->vectrix.resize(elem.size());
    for (unsigned i = 0; i < elem.size(); i++)
    {
        this->vectrix[i] = elem[i];
    }

    rows = elem.size();
    cols = 0;
}

template <typename T>
vx<T>::vx(unsigned _rows, unsigned _cols, const T &_initial)
{
    if (_cols == 0) // assigning a vector
    {
        unsigned sz = _rows;
        this->vectrix.resize(sz, _initial);
        rows = _rows;
        cols = _cols;
    }
    else
    {
        unsigned sz = _rows * _cols;
        this->vectrix.resize(sz, _initial);
        rows = _rows;
        cols = _cols;
    }
}

// duplicate constructor
template <typename T>
vx<T>::vx(const vx<T> &elem)
{
    vectrix = elem.vectrix;
    rows = elem.nrows();
    cols = elem.ncols();
}

// constructor from 2D array
template <typename T>
vx<T>::vx(const std::vector<std::vector<T>> &elem)
{

    unsigned r = elem.size();
    unsigned c = elem[0].size();
    unsigned len = r * c;
    this->vectrix.reserve(len);
    for (unsigned i = 0; i < r; i++)
    {
        for (unsigned j = 0; i < c; i++)
        {
            this->vectrix[GetIndex(cols, i, j)] = elem[i][j];
        }
    }

    rows = r;
    cols = c;
}

// virtual destructor
template <typename T>
vx<T>::~vx() {}

///// OPPORATORS/////:

//Matrix operations:

// assignment (=)
template <typename T>
vx<T> &vx<T>::operator=(const vx<T> &elem)
{
    if (&elem == this)
    {
        return *this;
    }
    else if (elem.ncols() == 0) // if it's a vector
    {
        unsigned newrows = elem.nrows();
        unsigned newcols = 0;
        vectrix.resize(newrows);
        for (unsigned i = 0; i < newrows; i++)
        {
            this->vectrix[i] = elem(i);
        }
        rows = newrows;
        cols = 0;
        return this->vectrix;
    }
    else if (elem.ncols() != 0) // it's a matrix
    {
        unsigned newrows = elem.nrows();
        unsigned newcols = elem.ncols();
        unsigned sz = newrows * newcols;

        vectrix.resize(sz);
        for (unsigned i = 0; i < newrows; i++)
        {
            for (unsigned j = 0; j < newcols; j++)
            {
                this->vectrix[GetIndex(newcols, i, j)] = elem(i, j);
            }
        }

        rows = newrows;
        cols = newcols;
        return this->vectrix;
    }
}

// matrix addition
template <typename T>
vx<T> vx<T>::operator+(const vx<T> &elem)
{
    vx result(rows, cols, 0.0);
    unsigned sz = rows * cols;
    unsigned sz2 = elem.nrows() * elem.ncols();
    if (sz != sz2) // if they're not equal size
    {
        std::cout << "* Operands could not be broadcasted together with their shapes *"
                  << "\n";
        exit(0);
    }

    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            result(i, j) = this->vectrix[GetIndex(cols, i, j)] + elem(i, j);
        }
    }
    return result;
}

// in-place matrix addition (cumulative addition; +=)
template <typename T>
vx<T> &vx<T>::operator+=(const vx<T> &elem)
{
    unsigned _rows = elem.nrows();
    unsigned _cols = elem.ncols();

    for (unsigned i = 0; i < _rows; i++)
    {
        for (unsigned j = 0; j < _cols; j++)
        {
            this->vectrix[GetIndex(cols, i, j)] += elem(i, j);
        }
    }
    return *this;
}

// matrix multiplication
template <typename T>
vx<T> matmul(const vx<T> &elem1, const vx<T> &elem2)
{
    unsigned m1_rows = elem1.nrows();
    unsigned m1_cols = elem1.ncols();
    unsigned m2_rows = elem2.nrows();
    unsigned m2_cols = elem2.ncols();

    if (m1_cols == m2_rows)
    {
        vx<T> result(m1_rows, m2_cols, 0.0);
        for (unsigned i = 0; i < m1_rows; i++)
        {
            for (unsigned j = 0; j < m2_cols; j++)
            {
                for (unsigned k = 0; k < m1_cols; k++)
                {
                    result(i, j) += elem1(i, k) * elem2(k, j);
                }
            }
        }
        return result;
    }
    else
    {
        std::cout << "* m1 cols != m2 rows *"
                  << "\n\n";
        exit(0);
    }
}

// element wise mutliplication - matrices
template <typename T>
vx<T> vx<T>::operator*(const vx<T> &elem)
{
    if (elem.ncols() == 0 && cols == 0) // both vecs
    {
        vx<T> result(cols, 0, 0);
        for (unsigned i = 0; i < this->vectrix.size(); i++)
        {
            result(i, 0) = this->vectrix[i] * elem(i, 0);
        }
        return result;
    }
    else if (elem.ncols() != 0 && cols != 0) // both matrices
    {
        unsigned m1_rows = rows;
        unsigned m1_cols = cols;
        unsigned m1_sz = m1_rows * m1_cols;
        unsigned m2_rows = elem.nrows();
        unsigned m2_cols = elem.ncols();
        unsigned m2_sz = m2_rows * m2_cols;
        vx<T> result(m1_rows, m1_cols, 0.0);
        if (m1_rows == m2_rows && m1_cols == m2_cols)
        {
            for (unsigned i = 0; i < m1_rows; i++)
            {
                for (unsigned j = 0; j < m1_cols; j++)
                {
                    result(i, j) = elem(i, j) * this->vectrix[GetIndex(cols, i, j)];
                }
            }
            return result;
        }
        else
        {
            std::cout << "* Elements have different shapes, cannot compute *"
                      << "\n\n";
            exit(0);
        }
    }
    else if (elem.ncols() == 0 && cols != 0) // elem is vec obj isnt
    {
        if (elem.nrows() == cols) // required for vec to mat multiplication
        {
            vx<T> result(rows, 0, 0);
            for (unsigned i = 0; i < rows; i++)
            {
                T rowDP = 0;
                for (unsigned j = 0; j < cols; j++)
                {
                    rowDP += this->vectrix[GetIndex(cols, i, j)] * elem(j, 0);
                }
                result(i, 0) = rowDP;
            }
            return result;
        }
        else
        {
            std::cout << "* Check inputs, cannot compute *"
                      << "\n";
        }
    }
    else if (elem.ncols() != 0 && cols == 0) // elem isnt vec obj is
    {
        if (elem.ncols() == rows) // required for vec to mat multiplication
        {
            vx<T> result(elem.nrows(), 0, 0);
            for (unsigned i = 0; i < elem.nrows(); i++)
            {
                T rowDP = 0;
                for (unsigned j = 0; j < elem.ncols(); j++)
                {
                    rowDP += this->vectrix[j] * elem(i, j);
                }
                result(i, 0) = rowDP;
            }
            return result;
        }
        else
        {
            std::cout << "* Check inputs, cannot compute *"
                      << "\n";
        }
    }
    else
    {
        std::cout << "* Cannot compute computation, check inputs *"
                  << "\n";
        exit(0);
    }
}

// matrix cumulative multiplication
template <typename T>
vx<T> &vx<T>::operator*=(const vx<T> &elem)
{
    vx result = (*this) * elem;
    *this = result;
    return *this;
}

// matrix subtraction
template <typename T>
vx<T> vx<T>::operator-(const vx<T> &elem)
{
    vx result(rows, cols, 0.0);

    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            result(i, j) = this->vectrix[GetIndex(cols, i, j)] - elem(i, j);
        }
    }
    return result;
}

template <typename T>
vx<T> &vx<T>::operator-=(const vx<T> &elem)
{
    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            this->vectrix[GetIndex(cols, i, j)] -= elem(i, j);
        }
    }
    return *this;
}

// SCALAR OPERATIONS TO MATRICES AND VECTORS

template <typename T>
vx<T> vx<T>::operator+(const T elem)
{
    vx result(rows, cols, 0.0);
    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            result(i, j) = this->vectrix[GetIndex(cols, i, j)] + elem;
        }
    }

    return result;
}

template <typename T>
vx<T> &vx<T>::operator+=(const T elem)
{
    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            this->vectrix[GetIndex(cols, i, j)] += elem;
        }
    }
    return *this;
}

template <typename T>
vx<T> vx<T>::operator-(const T elem)
{
    vx result(rows, cols, 0.0);

    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            result(i, j) = this->vectrix[GetIndex(cols, i, j)] - elem;
        }
    }

    return result;
}

template <typename T>
vx<T> &vx<T>::operator-=(const T elem)
{
    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            this->vectrix[GetIndex(cols, i, j)] -= elem;
        }
    }
    return *this;
}

template <typename T>
vx<T> vx<T>::operator*(const T elem)
{
    vx result(rows, cols, 0.0);

    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            result(i, j) = this->vectrix[GetIndex(cols, i, j)] * elem;
        }
    }

    return result;
}

template <typename T>
vx<T> &vx<T>::operator*=(const T elem)
{
    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            this->vectrix[GetIndex(cols, i, j)] *= elem;
        }
    }
    return *this;
}

template <typename T>
vx<T> vx<T>::operator/(const T elem)
{
    vx result(rows, cols, 0.0);

    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            result(i, j) = this->vectrix[GetIndex(cols, i, j)] / elem;
        }
    }

    return result;
}

template <typename T>
vx<T> &vx<T>::operator/=(const T elem)
{
    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            this->vectrix[GetIndex(cols, i, j)] /= elem;
        }
    }
    return *this;
}

template <typename T>
vx<T> vx<T>::operator^(const T elem)
{
    vx result(rows, cols, 0.0);

    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            result(i, j) = powf(this->vectrix[GetIndex(cols, i, j)], elem);
        }
    }

    return result;
}

// VECTOR TO VECTOR OPERATIONS

template <typename T>
std::vector<T> &add(std::vector<T> &elem1, std::vector<T> &elem2) // add two vectors
{
    unsigned v1_nitems = elem1.size();
    unsigned v2_nitems = elem2.size();

    if (v1_nitems > v2_nitems)
    {
        for (unsigned i = 0; i < v2_nitems; i++)
        {
            elem1[i] += elem2[i];
        }
        return elem1;
    }
    else if (v1_nitems < v2_nitems)
    {
        for (unsigned i = 0; i < v1_nitems; i++)
        {
            elem2[i] += elem1[i];
        }
        return elem2;
    }
    else
    {
        for (unsigned i = 0; i < v1_nitems; i++)
        {
            elem1[i] += elem2[i];
        }
        return elem1;
    }
}

template <typename T>
std::vector<T> &add(std::vector<T> &elem1, const T &elem2)
{
    for (unsigned i = 0; i < elem1.size(); i++)
    {
        elem1[i] += elem2;
    }
    return elem1;
}

template <typename T>
std::vector<T> &subtract(std::vector<T> &elem1, std::vector<T> &elem2)
{
    unsigned v1_nitems = elem1.size();
    unsigned v2_nitems = elem2.size();

    if (v1_nitems > v2_nitems)
    {
        for (unsigned i = 0; i < v2_nitems; i++)
        {
            elem1[i] -= elem2[i];
        }
        return elem1;
    }
    else if (v1_nitems < v2_nitems)
    {
        for (unsigned i = 0; i < v1_nitems; i++)
        {
            elem2[i] -= elem1[i];
        }
        return elem2;
    }
    else
    {
        for (unsigned i = 0; i < v1_nitems; i++)
        {
            elem1[i] -= elem2[i];
        }
        return elem1;
    }
}

template <typename T>
std::vector<T> &subtract(std::vector<T> &elem1, const T &elem2)
{
    for (unsigned i = 0; i < elem1.size(); i++)
    {
        elem1[i] -= elem2;
    }
    return elem1;
}

template <typename T>
std::vector<T> &scale(std::vector<T> &elem1, const T &elem2)
{
    for (unsigned i = 0; i < elem1.size(); i++)
    {
        elem1[i] *= elem2;
    }
    return elem1;
}

template <typename T>
float vx<T>::dot_product(const vx<T> &elem)
{ // both are vectors and their dims equal
    if (cols == 0 && elem.ncols() == 0 && rows == elem.nrows())
    {
        float dproduct = 0;
        for (unsigned i = 0; i < rows; i++)
        {
            dproduct += this->vectrix[i] * elem(i, 0);
        }
        return dproduct;
    }
    else
    {
        std::cout << "* Incompatable shapes: reshape the vectors to compute *"
                  << "\n\n";
        exit(0);
    }
}

template <typename T>
float vx<T>::cross_product(const vx<T> &elem)
{
    float v1_x = this->vector.size();
    float v1_y = this->vector[this->vector.size() - 1] - this->vector[0];
    float v2_x = elem.nrows();
    float v2_y = elem(elem.nrows() - 1) - elem(0);
    float dproduct = v1_x * v2_x + v1_y * v2_y;
    return dproduct;
}

template <typename T>
vx<T> expand(const std::vector<T> &elem)
{
    unsigned n = elem.size();
    vx<T> result(1, n, 0.0);
    for (unsigned i = 0; i < n; i++)
    {
        result(0, i) = elem[i];
    }
    return result;
}

template <typename T>
float angle_between(const std::vector<T> &elem1, const std::vector<T> &elem2) // measured in degrees
{
    float DP = dot_product(elem1, elem2);
    float mag1 = distance(elem1);
    float mag2 = distance(elem2);
    float mag_product = mag1 * mag2;
    float rad = deg_to_rad(DP / mag_product);
    float angle = acosf(rad);
    return angle;
}

template <typename T>
float distance(const std::vector<T> &elem)
{
    float startx = 1;
    float endx = elem.size();
    float starty = elem[0];
    float endy = elem[elem.size() - 1];
    float delta_x_sq = powf(endx - startx, 2);
    float delta_y_sq = powf(endy - starty, 2);
    return sqrtf(delta_x_sq + delta_y_sq);
}

template <typename T>
float vec_angle(const std::vector<T> &elem) // measured in degrees
{
    float slope1 = slope(elem);
    float direction = atanf(slope1);
    return direction;
}

template <typename T>
float slope(const std::vector<T> &elem)
{
    float startx = 1;
    float endx = elem.size();
    float starty = elem[0];
    float endy = elem[elem.size() - 1];
    float delta_y = endy - starty;
    float delta_x = endx - startx;
    float slope = delta_y / delta_x;
    return slope;
}

template <typename T>
std::vector<T> vx<T>::diag_vec()
{
    std::vector<T> result(rows, 0.0);
    for (unsigned i = 0; i < rows; i++)
    {
        result[i] = this->matrix[i][i];
    }
    return result;
}

template <typename T>
vx<T> vx<T>::transpose()
{
    unsigned old_rows = rows;
    unsigned old_cols = cols;
    unsigned new_rows = old_cols;
    unsigned new_cols = old_rows;
    vx result(new_rows, new_cols, 0.0);
    unsigned rows_curr = 0;
    unsigned cols_curr = 0;

    for (unsigned i = 0; i < old_rows; i++)
    {
        for (unsigned j = 0; j < old_cols; j++)
        {
            if (cols_curr == new_cols)
            {
                rows_curr++;
                cols_curr = 0;
            }
            result(rows_curr, cols_curr) = this->matrix[i][j];
            cols_curr++;
        }
    }
    return result;
}

// sum of matrix
template <typename T>
float vx<T>::sum()
{
    /*
    Returns: sum of matrix
    */
    float mx_sum = 0;
    for (int i = 0; i < this->matrix.size(); i++)
    {
        for (int j = 0; j < this->matrix[i].size(); j++)
            mx_sum += this->matrix[i][j];
    }
    return mx_sum;
}

// sum of vector
template <typename T>
float sum(const std::vector<T> &elem)
{
    float sum = 0.0;
    unsigned r = elem.size();
    for (unsigned i = 0; i < r; i++)
    {
        sum += elem[i];
    }
    return sum;
}

// sum of 2D vector
template <typename T>
float sum(const std::vector<std::vector<T>> &elem)
{
    float mx_sum = 0;
    for (int i = 0; i < elem.size(); i++)
    {
        for (int j = 0; j < elem[i].size(); j++)
            mx_sum += elem[i][j];
    }
    return mx_sum;
}

// mean of matrix
template <typename T>
float vx<T>::mean()
{
    /*
    Returns: mean of matrix
    */

    float sum = 0;
    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            sum += this->matrix[i][j];
        }
    }
    unsigned nitems = rows * cols;
    float mean = sum / nitems;
    return mean;
}

// mean of vector
template <typename T>
float mean(const std::vector<T> &elem)
{
    float sum1 = sum(elem);
    unsigned nitems = elem.size();
    float mean = sum1 / nitems;
    return mean;
}

// mean of 2D vector
template <typename T>
float mean(const std::vector<std::vector<T>> &elem)
{
    float sum = vx<T>::sum(elem);
    unsigned nitems = vx<T>::nitems(elem);
    float mean = sum / nitems;
    return mean;
}

// Standard deviation of a matrix
template <typename T>
float vx<T>::SD()
{
    float sq_var = 0; // sqaured variance term
    float sum = 0;
    unsigned nitems = rows * cols;
    for (unsigned n = 0; n < this->matrix.size(); n++)
    {
        for (unsigned z = 0; z < this->matrix[n].size(); z++)
        {
            sum += this->matrix[n][z];
        }
    }

    float mean = sum / nitems;
    for (unsigned i = 0; i < this->matrix.size(); i++)
    {
        for (unsigned j = 0; j < this->matrix[i].size(); j++)
        {
            sq_var += powf(this->matrix[i][j] - mean, 2);
        }
    }
    return sqrtf(sq_var / nitems);
}

// Standard deviation of a vector
template <typename T>
float SD(std::vector<T> &elem)
{
    float sq_var = 0; // sqaured variance term
    float mean1 = mean(elem);
    unsigned nitems1 = elem.size();
    for (unsigned i = 0; i < elem.size(); i++)
    {
        sq_var += powf(elem[i] - mean1, 2);
    }

    return sqrtf(sq_var / nitems1);
}

// Standard deviation of a 2D vector
template <typename T>
float SD(std::vector<std::vector<T>> &elem)
{
    float sq_var = 0; // sqaured variance term
    float mean1 = mean(elem);
    unsigned nitems1 = nitems(elem);
    for (unsigned i = 0; i < elem.size(); i++)
    {
        for (unsigned j = 0; j < elem[i].size(); j++)
        {
            sq_var += powf(elem[i][j] - mean1, 2);
        }
    }
    return sqrtf(sq_var / nitems1);
}

// flatten a matrix into a 1D vector
template <typename T>
std::vector<T> vx<T>::flatten()
{
    std::vector<T> new_vec;
    new_vec.reserve(rows * cols);

    for (unsigned i = 0; i < this->matrix.size(); i++)
    {
        for (unsigned j = 0; j < this->matrix[i].size(); j++)
        {
            new_vec.push_back(this->matrix[i][j]);
        }
    }
    return new_vec;
}

// max of a matrix
template <typename T>
T vx<T>::max()
{
    T item = this->matrix[0][0];
    for (unsigned i = 0; i < this->matrix.size(); i++)
    {
        for (unsigned j = 0; j < this->matrix[i].size(); j++)
        {
            if (this->matrix[i][j] > item)
            {
                item = this->matrix[i][j];
            }
        }
    }
    return item;
}

// max of a vector
template <typename T>
T max(const std::vector<T> &elem)
{
    T item = elem[0];
    for (unsigned i = 0; i < elem.size(); i++)
    {
        if (elem[i] > item)
        {
            item = elem[i];
        }
    }
    return item;
}

// max of a 2D vector
template <typename T>
T max(const std::vector<std::vector<T>> &elem)
{
    T item = elem[0][0];
    for (unsigned i = 0; i < elem.size(); i++)
    {
        for (unsigned j = 0; j < elem[i].size(); j++)
        {
            if (elem[i][j] > item)
            {
                item = elem[i][j];
            }
        }
    }
    return item;
}

// min of matrix
template <typename T> // generic variable type initialization
T vx<T>::min()
{
    /*
        Gets the minimum value of matrix
        Parameter "matrix": 2D matrix
        Returns: min value of matrix
        */
    T item = this->matrix[0][0];
    for (int i = 0; i < this->matrix.size(); i++)
    {
        for (int j = 0; j < this->matrix[i].size(); j++)
            if (this->matrix[i][j] < item)
            {
                item = this->matrix[i][j];
            }
    }
    return item;
}

// min of vector
template <typename T>
T min(const std::vector<T> &elem)
{
    T item = elem[0];
    for (int i = 0; i < elem.size(); i++)
    {
        if (elem[i] < item)
        {
            item = elem[i];
        }
    }
    return item;
}

// min of 2D vector
template <typename T>
T min(const std::vector<std::vector<T>> &elem)
{
    T item = elem[0][0];
    for (int i = 0; i < elem.size(); i++)
    {
        for (int j = 0; j < elem[i].size(); j++)
            if (elem[i][j] < item)
            {
                item = elem[i][j];
            }
    }
    return item;
}

// abs values of a matrix
template <typename T>
vx<T> &vx<T>::absolute()
{
    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            this->matrix[i][j] = sqrtf(powf(this->matrix[i][j], 2));
        }
    }
    return *this;
}

// abs value of vector
template <typename T>
std::vector<T> &absolute(std::vector<T> &elem)
{
    unsigned r = elem.size();

    for (unsigned i = 0; i < r; i++)
    {
        elem[i] = sqrtf(powf(elem[i], 2));
    }
    return elem;
}

// abs of a 2D vector
template <typename T>
std::vector<std::vector<T>> &absolute(std::vector<std::vector<T>> &elem)
{
    unsigned r = elem.size();
    unsigned c = elem[0].size();

    for (unsigned i = 0; i < r; i++)
    {
        for (unsigned j = 0; j < c; j++)
        {
            elem[i][j] = sqrtf(powf(elem[i][j], 2));
        }
    }
    return elem;
}

template <typename T>
T &absolute(T value)
{
    return sqrtf(powf(value, 2));
}

// balances two matrices to be of the same order (dims)
template <typename T>
std::vector<vx<T>> vx<T>::balance(const vx<T> &elem1, const vx<T> &elem2, T value)
{
    unsigned m1_rows = elem1.ncols();
    unsigned m2_rows = elem2.nrows();
    unsigned m1_cols = elem1.ncols();
    unsigned m2_cols = elem2.ncols();
    std::vector<unsigned> rvec = {m1_rows, m2_rows};
    std::vector<unsigned> cvec = {m1_cols, m2_cols};
    unsigned r = max(rvec);
    unsigned c = max(cvec);
    vx new_mx(r, 1, value);
    vx temp_vec(1, value);

    if (m1_rows > m2_rows)
    {
        // filling mx with shorter rows so that they're equal row number
        for (unsigned n = 0; n < m1_rows - m2_rows; n++)
        {
            elem2.insert(elem2.end(), temp_vec);
        }

        for (unsigned i = 0; i < r; i++) // loop through rows
        {
            if (elem2[i].size() > elem1[i].size())
            {
                for (unsigned n = 0; n = elem2[i].size() - elem1[i].size(); n++) // appending 0 to shorter mx cols wise (m1 in this case) for i times
                {
                    elem1[i].push_back(value);
                }
            }
            else if (elem2[i].size() < elem1[i].size())
            {
                for (unsigned n = 0; n = elem1[i].size() - elem2[i].size(); n++) // appending 0 to shorter mx cols wise (m1 in this case) for i times
                {
                    elem2[i].push_back(value);
                }
            }
        }
    }
    else if (m1_rows < m2_rows)
    {
        // filling mx with shorter rows so that they're equal row number
        for (unsigned n = 0; n < m2_rows - m1_rows; n++)
        {
            elem1.insert(elem1.end(), temp_vec);
        }
        // working with m2
        for (unsigned i = 0; i < r; i++) // loop through rows of m2
        {
            if (elem2[i].size() > elem1[i].size())
            {
                for (unsigned n = 0; n = elem2[i].size() - elem1[i].size(); n++) // appending 0 to shorter mx cols wise (m1 in this case) for i times
                {
                    elem1[i].push_back(value);
                }
            }
            else if (elem2[i].size() < elem1[i].size())
            {
                for (unsigned n = 0; n = elem1[i].size() - elem2[i].size(); n++) // appending 0 to shorter mx cols wise (m1 in this case) for i times
                {
                    elem2[i].push_back(value);
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < r; i++) // loop through rows of m2
        {
            if (elem2[i].size() > elem1[i].size())
            {
                for (int n = 0; n = elem2[i].size() - elem1[i].size(); n++) // appending 0 to shorter mx cols wise (m1 in this case) for i times
                {
                    elem1[i].push_back(value);
                }
            }
            else if (elem2[i].size() < elem1[i].size())
            {
                for (int n = 0; n = elem1[i].size() - elem2[i].size(); n++) // appending 0 to shorter mx cols wise (m1 in this case) for i times
                {
                    elem2[i].push_back(value);
                }
            }
        }
    }
    std::vector<vx<T>> balanced = {elem1, elem2};
    return balanced;
}

// UTLITY FUNCTIONS:

template <typename T>
std::vector<unsigned> vx<T>::dims(bool print) const
{
    unsigned H = rows;
    unsigned L = cols;
    if (print == true)
    {
        std::cout << "Dimensions: [" << H << ", " << L << "] \n";
    }
    std::vector<unsigned> dims = {H, L};
    return dims;
}

template <typename T>
std::vector<unsigned> dims(const std::vector<std::vector<T>> &elem, bool print)
{
    unsigned H = elem.size();
    unsigned L = elem[0].size();
    if (print == true)
    {
        std::cout << "Dimensions: [" << H << ", " << L << "] \n";
    }
    std::vector<unsigned> dims = {H, L};
    return dims;
}

template <typename T>
unsigned vx<T>::nitems() const
{
    if (cols == 0)
    {
        return rows;
    }
    else
    {
        return rows * cols;
    }
}

template <typename T>
unsigned nitems(const std::vector<std::vector<T>> &elem)
{
    unsigned r = nrows(elem);
    unsigned c = ncols(elem);
    return r * c;
}

template <typename T>
unsigned vx<T>::nrows() const
{
    return this->rows;
}

template <typename T>
unsigned nrows(const std::vector<std::vector<T>> &elem)
{
    unsigned r = elem.size();
    return r;
}

template <typename T>
unsigned vx<T>::ncols() const
{
    return this->cols;
}

template <typename T>
unsigned ncols(const std::vector<std::vector<T>> &elem)
{
    unsigned c = elem[0].size();
    return c;
}

template <typename T>
std::vector<T> range1D(const T stop, const T start, const T increment)
{
    std::vector<T> vec;
    for (T i = start; i <= stop; i += increment)
    {
        vec.push_back(i);
    }
    return vec;
}

template <typename T>
vx<T> range2D(const T stop, const T start, const T increment, const unsigned numrows)
{
    std::vector<T> vec = range1D(stop, start, increment);
    std::vector<std::vector<T>> mx(numrows, vec);
    vx<T> result(mx);
    return result;
}

template <typename T>
vx<T> vx<T>::reshape(const unsigned rows, const unsigned cols)
{
    unsigned old_rows = this->rows;
    unsigned old_cols = this->cols;
    unsigned curr_rows = 0;
    unsigned curr_cols = 0;
    vx result(rows, cols, 0.0);

    for (int i = 0; i < old_rows; i++)
    {
        for (int j = 0; j < old_cols; j++)
        {
            if (curr_cols == cols)
            {
                curr_rows++;
                curr_cols = 0;
            }
            result(curr_rows, curr_cols) = this->matrix[i][j];
            curr_cols++;
        }
    }
    return result;
}

template <typename T>
float rad_to_deg(T value)
{
    const float one80_over_pi = 180 / pi;
    return value * one80_over_pi;
}

template <typename T>
float deg_to_rad(T value)
{
    const float pi_over_180 = pi / 180;
    return value * pi_over_180;
}

template <typename T>
bool is_ragged(const std::vector<std::vector<T>> &elem)
{
    unsigned num_colitems = elem[0].size();
    for (unsigned i = 0; i < elem.size(); i++)
    {
        if (elem[i].size() != num_colitems)
        {
            return true;
        }
    }
    return false;
}

template <typename T>
bool vx<T>::is_ragged()
{
    unsigned num_colitems = this->matrix[0].size();
    for (unsigned i = 0; i < this->matrix.size(); i++)
    {
        if (this->matrix[i].size() != num_colitems)
        {
            return true;
        }
    }
    return false;
}

template <typename T>
unsigned digit_count(const T &num)
{
    return unsigned(log10f(num) + 1);
}

template <typename T>
unsigned vx<T>::digit_count() const
{
    unsigned digits = 0;
    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            unsigned dig = unsigned(log10f(this->matrix[i][j]) + 1);
            if (dig > digits)
            {
                digits = dig;
            }
        }
    }
    return digits;
}

template <typename T>
unsigned digit_count(const std::vector<T> &elem)
{
    unsigned digits = 0;
    for (unsigned i = 0; i < elem.size(); i++)
    {
        unsigned dig = unsigned(log10f(elem[i]) + 1);
        if (dig > digits)
        {
            digits = dig;
        }
    }
    return digits;
}

template <typename T>
void vx<T>::print() const
{
    if (cols == 0)
    {
        std::cout << "[" << this->vectrix.size() << " Member Vector]:"
                  << "\n\n";
        std::cout << "[  ";
        for (unsigned i = 0; i < this->vectrix.size(); i++)
        {
            std::cout << this->vectrix[i] << "  ";
        }
        std::cout << "]";
        std::cout << "\n\n";
    }
    else
    {
        std::cout << "[" << rows << " x " << cols << "] "
                  << "Matrix: "
                  << "\n\n";
        for (unsigned i = 0; i < rows; i++)
        {
            std::cout << "[  ";
            for (unsigned j = 0; j < cols; j++)
            {
                std::cout << this->vectrix[GetIndex(cols, i, j)] << "  ";
            }
            std::cout << "]"
                      << "\n";
        }
        std::cout << "\n\n";
    }
}

template <typename T>
void print(const T &elem)
{
    std::cout << elem << "\n\n";
}

template <typename T>
void print(const std::vector<T> &elem)
{
    std::cout << "[" << elem.size() << " Member Vector]:"
              << "\n\n";
    std::cout << "[  ";
    for (unsigned i = 0; i < elem.size(); i++)
    {
        std::cout << elem[i] << "  ";
    }
    std::cout << "]";
    std::cout << "\n\n";
}

template <typename T>
void print(const std::vector<std::vector<T>> &elem)
{
    std::cout << "[" << elem.size() << " x " << elem[0].size() << "] "
              << "2D Vector: "
              << "\n\n";
    for (unsigned i = 0; i < elem.size(); i++)
    {
        std::cout << "[  ";
        for (unsigned j = 0; j < elem[0].size(); j++)
        {
            std::cout << elem[i][j] << "  ";
        }
        std::cout << "]"
                  << "\n";
    }
    std::cout << "\n\n";
}

template <typename T>
std::vector<T> &sort(std::vector<T> &elem)
{
    for (int k = 1; k < elem.size(); k++)
    {
        T t = elem[k];
        int j = k - 1;

        while (j >= 0 && t <= elem[j])
        {
            elem[j + 1] = elem[j];
            j = j - 1;
        }
        elem[j + 1] = t;
    }
    return elem;
}

template <typename T>
float median(std::vector<T> &elem)
{
    unsigned n = elem.size() + 1;
    unsigned loc = n / 2 - 1;
    std::vector<T> sorted = sort(elem);

    if (elem.size() % 2 == 0)
    {
        return (sorted[elem.size() / 2 - 1] + sorted[elem.size() / 2]) / 2;
    }

    return sorted[loc];
}

template <typename T>
float vx<T>::median() const
{
    unsigned n = rows * cols + 1;
    unsigned loc = n / 2 - 1;
    std::vector<T> new_vec;
    new_vec.reserve(rows * cols);
    for (unsigned i = 0; i < this->matrix.size(); i++)
    {
        for (unsigned j = 0; j < this->matrix[i].size(); j++)
        {
            new_vec.push_back(this->matrix[i][j]);
        }
    }
    std::vector<T> sorted = sort(new_vec);
    if (sorted.size() % 2 == 0)
    {
        return (sorted[new_vec.size() / 2 - 2] + sorted[new_vec.size() / 2]) / 2;
    }

    return new_vec[loc];
}

template <typename T>
T &vx<T>::operator()(const unsigned &row, const unsigned &col)
{
    if (cols == 0 || col == 0)
    {
        return this->vectrix[row];
    }
    else
    {
        return this->vectrix[GetIndex(cols, row, col)];
    }
}

template <typename T>
const T &vx<T>::operator()(const unsigned &row, const unsigned &col) const
{
    if (cols == 0 || col == 0)
    {
        return this->vectrix[row];
    }
    else
    {
        return this->vectrix[GetIndex(cols, row, col)];
    }
}

std::vector<double> getRandom_decimal(const unsigned n)
{
    std::vector<double> vec;
    vec.reserve(n);
    srand((unsigned)time(NULL));

    for (int i = 0; i < n; i++)
    {
        double thing = (double)rand() / RAND_MAX;
        vec.push_back(thing);
    }
    return vec;
}

#endif