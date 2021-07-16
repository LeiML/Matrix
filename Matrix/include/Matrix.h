//
// Created by LeiLei on 2021/6/26.
//

#ifndef MATRIX_MATRIX_H
#define MATRIX_MATRIX_H
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <random>

class Size{
public:
    Size(const int& row, const int& col);
    friend std::ostream &operator<< (std::ostream& os, const Size& size);
    bool operator== (const Size& size)const;
    bool operator!= (const Size& size)const;
    int row, col;
};

template<typename Tp>
class MatrixCommand{
public:
    MatrixCommand(const int& row, const int& col) : row(row), col(col) {}
    MatrixCommand<Tp>& operator, (const Tp& num);
    std::vector<Tp> vec{};
    int row, col;
};

template<typename Tp>
MatrixCommand<Tp> &MatrixCommand<Tp>::operator,(const Tp &num) {
    this->vec.push_back(num);
    return *this;
}

template<typename Tp>
class Matrix_{
public:
    Matrix_(const int& row, const int& col) : row(row), col(col) {}
    int row, col;
};

template <typename Tp, typename Tp2>inline
MatrixCommand<Tp> operator<<(const Matrix_<Tp>& mat, const Tp2& num) {
    MatrixCommand<Tp>matrixCommand = MatrixCommand<Tp>(mat.row, mat.col);
    return (matrixCommand, (Tp)num);
}

class Matrix{
public:
    /* construct function */
    Matrix() = default;
    Matrix(const int&row, const int& col);
    explicit Matrix(const Size& size);
    template<class Tp>Matrix(const int&row, const int& col, const std::vector<Tp>&vec);
    template<class Tp>Matrix(const Size& size, const std::vector<Tp>&vec);
    Matrix(const Matrix& mat);
    template<class Tp> Matrix(const MatrixCommand<Tp>& command);
    virtual ~Matrix();
    /* static function */
    template<class Tp>static Matrix fill(const int& row, const int& col, const Tp& num);
    template<class Tp>static Matrix fill(const Size& size, const Tp& num);
    static Matrix ones(const int& row, const int& col) ;
    static Matrix ones(const Size& size);
    static Matrix zeros(const int& row, const int& col);
    static Matrix zeros(const Size& size);
    static Matrix eye(const int& row);
    static Matrix normal(const int& row, const int& col);
    /* Operator overloading */
    friend std::ostream& operator<< (std::ostream& os, const Matrix& mat);
    Matrix operator* (const Matrix& mat) const;
    template<typename Tp>Matrix operator* (const Tp& num) const;
    Matrix operator+ (const Matrix& mat) const;
    template<typename Tp>Matrix operator+ (const Tp& num) const;
    Matrix operator- (const Matrix& mat) const;
    Matrix operator- () const;
    template<typename Tp>Matrix operator- (const Tp& num) const;
    template<typename Tp>Matrix operator/ (const Tp& num) const;
    template<typename Tp>Matrix operator^(const Tp& num) const;
    Matrix& operator= (const Matrix& mat);
    bool operator== (const Matrix& mat) const;
    bool operator!= (const Matrix& mat) const;
    /* function */
    [[nodiscard]] Size size() const;
    [[nodiscard]] double &at(const int& r, const int& c) const;
    [[nodiscard]] Matrix rows(const int& r) const;
    [[nodiscard]] Matrix cols(const int& c) const;
    [[nodiscard]] Matrix transpose() const;
    [[nodiscard]] Matrix T() const;
    [[nodiscard]] double sum() const;
    [[nodiscard]] Matrix dot(const Matrix& mat) const;
    [[nodiscard]] bool empty() const;
    template<class Tp>std::vector<Tp>toVector() const;
    [[nodiscard]] Matrix reshape(const int& r, const int& col) const;
public:
    int row{}, col{};
private:
    double *ptr{};
    void allocation();
    void release();
};

template<class Tp>
Matrix::Matrix(const int &row, const int &col, const std::vector<Tp> &vec) : row(row), col(col) {
    this->allocation();
    for (int i = 0; i < this->row; ++i) {
        for (int j = 0; j < this->col; ++j) {
            this->at(i, j) = (i * this->col + j >= vec.size() ? 0:double(vec.at(i * this->col + j)));
        }
    }
}

template<class Tp>
Matrix::Matrix(const Size &size, const std::vector<Tp> &vec) : row(size.row), col(size.col){
    new (this)Matrix(this->row, this->col, vec);
}

template<class Tp>
Matrix::Matrix(const MatrixCommand<Tp> &command) : row(command.row), col(command.col) {
    new (this)Matrix(this->row, this->col, command.vec);
}

template<class Tp>
Matrix Matrix::fill(const int &row, const int &col, const Tp &num) {
    std::vector<Tp> vec(row * col, num);
    return Matrix(row, col, vec);
}

template<class Tp>
Matrix Matrix::fill(const Size &size, const Tp &num) {
    return Matrix::fill(size.row, size.col, num);
}

template<typename Tp>
Matrix Matrix::operator*(const Tp &num) const {
    return (*this) * fill(this->row, this->col, num);
}

template<typename Tp>
Matrix Matrix::operator+(const Tp &num) const {
    return (*this) + fill(this->row, this->col, num);
}

template<typename Tp>
Matrix Matrix::operator-(const Tp &num) const {
    return (*this) - fill(this->row, this->col, num);
}

template<typename Tp>
Matrix Matrix::operator/(const Tp &num) const {
    assert(num != 0);
    return (*this) * fill(this->row, this->col, 1/(double)num);
}
template<typename Tp>
Matrix Matrix::operator^ (const Tp &num) const {
    assert(!this->empty());
    return pow(*this, num);
}

template<class Tp>
std::vector<Tp> Matrix::toVector() const {
    std::vector<Tp>vec;
    for (int i = 0; i < this->row; ++i)
        for (int j = 0; j < this->col; ++j)
            vec.push_back((Tp)this->at(i, j));
    return vec;
}

inline
Matrix operator+ (const double& num, const Matrix& mat){
    return mat + num;
}

inline
Matrix operator- (const double& num, const Matrix& mat){
    return num + (-mat);
}
inline
Matrix operator* (const double& num, const Matrix& mat){
    return mat * num;
}

inline
Matrix pow(const Matrix& mat, const double& num){
    Matrix result = Matrix(mat.row, mat.col);
    for (int i = 0; i < mat.row; ++i)
        for (int j = 0; j < mat.col; ++j)
            result.at(i, j) = std::pow(mat.at(i, j), num);
    return result;
}

inline
Matrix exp(const Matrix& mat){
    Matrix result = Matrix(mat.row, mat.col);
    for (int i = 0; i < mat.row; ++i)
        for (int j = 0; j < mat.col; ++j)
            result.at(i, j) = std::exp(mat.at(i, j));
    return result;
}

inline
Matrix log(const Matrix& mat){
    Matrix result = Matrix(mat.row, mat.col);
    for (int i = 0; i < mat.row; ++i)
        for (int j = 0; j < mat.col; ++j)
            result.at(i, j) = std::log(mat.at(i, j));
    return result;
}

inline
double min(const Matrix& mat){
    assert(!mat.empty());
    double num = mat.at(0, 0);
    for (int i = 0; i < mat.row; ++i)
        for (int j = 0; j < mat.col; ++j)
            num = mat.at(i, j) < num ? mat.at(i, j) : num;
    return num;
}

inline
double max(const Matrix& mat){
    assert(!mat.empty());
    double num = mat.at(0, 0);
    for (int i = 0; i < mat.row; ++i)
        for (int j = 0; j < mat.col; ++j)
            num = mat.at(i, j) > num ? mat.at(i, j) : num;
    return num;
}

double mean(const Matrix& mat);

Matrix row_stack(const Matrix& mat, const int& r);

Matrix col_stack(const Matrix& mat, const int& c);

#endif //MATRIX_MATRIX_H
