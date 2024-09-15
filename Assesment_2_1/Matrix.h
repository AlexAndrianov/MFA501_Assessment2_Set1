#pragma once
#include <stdexcept>
#include <vector>
#include <math.h>
#include <iostream>

namespace math {

class Vector {
public:
    Vector() = default;
    Vector(std::vector<double> &&v) : _v(std::move(v)) {}

    Vector(size_t sz)
    {
        _v.resize(sz, 0.);
    }

    double length() const {
        double res = 0;

        for(auto i = 0u; i < _v.size(); i++)
            res += _v[i] * _v[i];

        return sqrtf(res);
    }

    Vector normV() const
    {
        auto res = *this;
        const auto len = length();

        for(auto i = 0u; i < res.size(); i++)
            res[i] /= len;

        return res;
    };

    double dot(const Vector &other) const
    {
        if(_v.size() != other._v.size())
            throw std::runtime_error("Size of vectors is not equial");

        auto res = 0.;
        for(auto i = 0u; i < size(); i++)
            res += _v[i] * other._v[i];

        return res;
    };

    double &operator[](size_t i) {
        return _v[i];
    }

    const double &operator[](size_t i) const {
        return _v[i];
    }

    Vector operator * (const double k) const
    {
        Vector res = *this;
        for(auto i = 0u; i < size(); i++)
            res[i] *= k;

        return res;
    }

    Vector &operator -= (const Vector &other)
    {
        if(_v.size() != other._v.size())
            throw std::runtime_error("Size of vectors is not equial");

        for(auto i = 0u; i < size(); i++)
            _v[i] -= other._v[i];

        return *this;
    }

    friend std::ostream& operator<<(std::ostream& out, const Vector& vc)
    {
        for(auto i = 0u; i < vc.size(); i++)
            out << vc._v[i] << std::endl;

        out << std::endl;
        return out;
    }

    size_t size() const { return _v.size(); }

    std::vector<double> _v;
};

class Matrix {
public:
    Matrix() = default;

    Matrix(size_t row, size_t col)
    {
        if(!row || !col)
            throw std::runtime_error("Empty matrix");

        _v.resize(row, Vector(col));
    }

    void addRow(Vector &&row)
    {
        if(!_v.empty() && _v[0].size() != row.size())
            throw std::runtime_error("Numbers of columns is defferent");

        _v.emplace_back(std::move(row));
    }

    bool isEmpty() const
    {
        return _v.empty();
    }

    bool isSquare() const
    {
        if(_v.empty())
            return true;

        return _v.size() == _v[0].size();
    }

    double &operator()(size_t r, size_t c) {
        return _v[r][c];
    }

    const double &operator()(size_t r, size_t c) const {
        return _v[r][c];
    }

    friend std::ostream& operator<<(std::ostream& out, const Matrix& mtx)
    {
        for(auto i = 0u; i < mtx._v.size(); i++)
            out << mtx._v[i];

        out << std::endl;
        return out;
    }

    Matrix operator * (const Matrix &mtrB) const {
        const auto &mtrA = *this;

        if(mtrA.isEmpty() || mtrB.isEmpty())
        {
            if(mtrB.isEmpty() != mtrA.isEmpty())
                throw std::runtime_error("Matrices have different sizes");

            return Matrix{};
        }

        const auto dotLenght = mtrA._v[0].size();

        if(dotLenght != mtrB._v.size())
            throw std::runtime_error("Matrices are not compatible");

        Matrix res{mtrA._v.size(), mtrB._v[0].size()};

        for(auto i = 0u; i < mtrA._v.size(); i++)
            for(auto j = 0u; j < mtrB._v[0].size(); j++)
                for(auto k = 0u; k < dotLenght; k++)
                    res(i,j) += mtrA._v[i][k] * mtrB._v[k][j];

        return res;
    }

    std::pair<Matrix, Matrix> decompositionQR() const
    {
        if(isEmpty())
            return {};

        if(!isSquare())
            throw std::runtime_error("The matrix should be square");

        Matrix Q{_v.size(), _v.size()};
        Matrix R{_v.size(), _v.size()};

        auto vectorByColumn = [](const auto &m, size_t col){
            if(m.isEmpty())
                throw std::runtime_error("The matrix is empty");

            if(col >= m._v[0].size())
                throw std::runtime_error("The column number exceeds matix size");

            Vector res(m._v.size());
            for(auto i = 0u; i < m._v.size(); i++)
                res[i] = m(i,col);

            return res;
        };

        Vector a_perp = vectorByColumn(*this, 0);

        for(auto i = 0u; i < _v.size(); i++)
        {
            const auto norm_perp = a_perp.length();
            const auto normV_perp = a_perp.normV();

            R(i, i) = norm_perp;

            for(auto j = 0u; j < normV_perp.size(); j++)
                Q(j, i) = normV_perp[j];

            const auto iNext = i + 1;
            if(iNext == _v.size())
                break;

            const auto ai = vectorByColumn(*this, iNext);
            auto a_perp_next = ai;

            for(auto j = 0u; j < iNext; j++)
            {
                const auto qj = vectorByColumn(Q, j);
                const auto k = ai.dot(qj);
                R(j, iNext) = k;
                a_perp_next -= qj * k;
            }

            a_perp = a_perp_next;
        }

        return {std::move(Q), std::move(R)};
    }

    bool isUpperTriangle(const double tol) const
    {
        if(!isSquare())
            throw std::runtime_error("The matrix should be square");

        double summ = 0.;
        for(auto j = 0; j < _v.size(); j++)
            for(auto i = j + 1; i < _v.size(); i++)
                summ += fabs(_v[i][j]);

        return summ < tol;
    }

    Vector eigenValuesByQRAlgorithm(const double tol = 1e-6, const size_t maxIterNumber = 1000) const
    {
        if(!isSquare())
            throw std::runtime_error("The matrix should be square");

        if(isEmpty())
            return {};

        Matrix a = *this;

        for(auto i = 0; i < maxIterNumber; i++)
        {
            const auto[q, r] = a.decompositionQR();

            std::cout << "Iteration n: " << i << std::endl;
            std::cout << "Q matrix: " << std::endl;
            std::cout << q;

            std::cout << "R matrix: " << std::endl;
            std::cout << r;

            a = r * q;

            std::cout << "Next A matrix: " << std::endl;
            std::cout << a;

            if(a.isUpperTriangle(tol))
                break;
        }

        Vector res(a._v.size());
        for(auto i = 0; i < a._v.size(); i++)
            res[i] = a(i, i);

        return res;
    }

    std::vector<Vector> _v;
};

}
