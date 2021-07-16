//
// Created by LeiLei on 2021/6/26.
//

#include "Matrix.h"

Size::Size(const int &row, const int &col) : row(row), col(col){}

std::ostream &operator<<(std::ostream &os, const Size &size) {
    os << "Size([" << size.row << ", " << size.col << "])";
    return os;
}

bool Size::operator==(const Size &size) const {
    return (this->row == size.row and this->col == size.col);
}

bool Size::operator!=(const Size &size) const {
    return (this->row != size.row or this->col != size.col);
}


Matrix::Matrix(const int &row, const int &col) : row(row), col(col){
    this->allocation();
}

Matrix::Matrix(const Size &size) : row(size.row), col(size.col){
    new (this)Matrix(this->row, this->col);
}

Matrix::Matrix(const Matrix& mat) : row(mat.row), col(mat.col) {
    this->allocation();
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            this->at(i, j) = mat.at(i, j);
        }
    }
}

Matrix::~Matrix() {
    this->release();
}

std::ostream &operator<<(std::ostream &os, const Matrix& mat) {
    assert(mat.row > 0 or mat.col > 0);
    os << "Matrix([";
    for (int i = 0; i < mat.row; ++i) {
        for (int j = 0; j < mat.col; ++j) {
            os << mat.at(i, j) << (j == mat.col - 1 ? "":", ");
        }
        os << (i == mat.row - 1 ? "])" : "\n        ");
    }
    return os;
}

void Matrix::allocation() {
    this->ptr = new double[this->row * this->col + 1];
}

void Matrix::release() {
    delete[] this->ptr;
}

double &Matrix::at(const int &r, const int &c) const {
    assert(this->row >= r and this->col >= c);
    return this->ptr[r * this->col + c];
}

Size Matrix::size() const {
    return Size(this->row, this->col);
}

Matrix Matrix::ones(const int &row, const int &col) {
    return fill(row, col, 1);
}

Matrix Matrix::ones(const Size &size) {
    return fill(size, 1);
}

Matrix Matrix::zeros(const int &row, const int &col) {
    return fill(row, col, 0);
}

Matrix Matrix::zeros(const Size &size) {
    return fill(size, 0);
}

Matrix Matrix::eye(const int &row) {
    Matrix mat = zeros(row, row);
    for (int i = 0; i < row; ++i)
        mat.at(i ,i) = 1;
    return mat;
}

Matrix Matrix::normal(const int& row, const int& col) {
    std::random_device engine;
    std::normal_distribution<double> distribution(0, 1);
    std::vector<double>vec;
    vec.reserve(row * col);
    for (int i = 0; i < row * col; ++i) {
        vec.push_back(distribution(engine));
    }
    return Matrix(row, col, vec);
}

Matrix Matrix::operator*(const Matrix &mat) const {
    assert(this->size() == mat.size());
    Matrix result = Matrix(this->size());
    for (int i = 0; i < this->row; ++i) {
        for (int j = 0; j < this->col; ++j) {
            result.at(i, j) = this->at(i, j) * mat.at(i, j);
        }
    }
    return result;
}

Matrix Matrix::operator+(const Matrix &mat) const {
    assert(this->size() == mat.size());
    Matrix result = Matrix(this->size());
    for (int i = 0; i < this->row; i++)
        for (int j = 0; j < this->col; ++j)
            result.at(i, j) = this->at(i, j) + mat.at(i, j);
    return result;
}

Matrix Matrix::operator-(const Matrix &mat) const {
    assert(this->size() == mat.size());
    Matrix result = Matrix(this->size());
    for (int i = 0; i < this->row; i++)
        for (int j = 0; j < this->col; ++j)
            result.at(i, j) = this->at(i, j) - mat.at(i, j);
    return result;
}

Matrix Matrix::operator-() const {
    return (*this) * fill(this->row, this->col, -1);
}

Matrix Matrix::rows(const int& r) const{
    Matrix mat = Matrix(1, this->col);
    for (int i = 0; i < this->col; ++i)
        mat.at(0, i) = this->at(r, i);
    return mat;
}

Matrix Matrix::cols(const int &c) const {
    Matrix mat = Matrix(this->row, 1);
    for (int i = 0; i < this->row; ++i)
        mat.at(i, 0) = this->at(i, c);
    return mat;
}

Matrix Matrix::transpose() const {
    Matrix mat = Matrix(this->col, this->row);
    for (int i = 0; i < this->row; ++i) {
        for (int j = 0; j < this->col; ++j) {
            mat.at(j, i) = this->at(i, j);
        }
    }
    return mat;
}

Matrix Matrix::dot(const Matrix& mat) const {
    assert(this->col == mat.row);
    Matrix result = Matrix(this->row, mat.col);
    for (int i = 0; i < this->row; ++i)
        for (int j = 0; j < mat.col; ++j)
            result.at(i, j) = (this->rows(i) * mat.cols(j).transpose()).sum();
    return result;
}

double Matrix::sum() const {
    double num = 0.;
    for (int i = 0; i < this->row; ++i)
        for (int j = 0; j < this->col; ++j)
            num += this->at(i, j);
    return num;
}

Matrix Matrix::T() const {
    return this->transpose();
}

bool Matrix::empty() const {
    return (this->row <=0 or this->col <= 0);
}

Matrix& Matrix::operator=(const Matrix &mat) {
    this->row = mat.row;
    this->col = mat.col;
    this->allocation();
    for (int i = 0; i < this->row; ++i)
        for (int j = 0; j < this->col; ++j)
            this->at(i, j) = mat.at(i, j);
    return *this;
}

Matrix Matrix::reshape(const int& r, const int& c) const {
    int m_r = r, m_c = c;
    if (r < 0) {
        assert((this->row * this->col) / (double) c - int((this->row * this->col) / (double) c) == 0);
        m_r = int((this->row * this->col) / (double) c);
    }
    if (c < 0) {
        assert((this->row * this->col) / (double) r - int((this->row * this->col) / (double) r) == 0);
        m_c = int((this->row * this->col) / (double) r);
    }
    assert(m_r * m_c == this->row * this->col);
    return Matrix(m_r, m_c, this->toVector<double>());
}

bool Matrix::operator==(const Matrix &mat) const {
    assert(!this->empty() and !mat.empty());
    if (this->size() != mat.size())
        return false;
    return this->toVector<double>() == mat.toVector<double>();
}

bool Matrix::operator!=(const Matrix &mat) const {
    return !(*this == mat);
}


double mean(const Matrix& mat){
    return mat.sum() / (mat.row * mat.col);
}

Matrix row_stack(const Matrix& mat, const int& r) {
    std::vector<double>vec;
    for (int i = 0; i < r; ++i) {
        std::vector<double>temp = mat.toVector<double>();
        vec.insert(vec.end(), temp.begin(), temp.end());
    }
    return Matrix(mat.row * r, mat.col, vec);
}

Matrix col_stack(const Matrix& mat, const int& c){
    std::vector<double>vec;
    for (int j = 0; j < mat.row; ++j) {
        for (int i = 0; i < c; ++i) {
            std::vector<double>temp = mat.rows(j).toVector<double>();
            vec.insert(vec.end(), temp.begin(), temp.end());
        }
    }
    return Matrix(mat.row, mat.col * c, vec);
}
