/* 
 * File:   MathMatrix.h
 * Author: tombieli
 *
 * Created on 6 November 2015, 15:22
 */

#ifndef MATHMATRIX_H
#define	MATHMATRIX_H

#include "MathVector.h"
#include <stdexcept>

namespace ksg
{

namespace utils
{
    
    namespace math
    {
        class MatrixInversionException : public std::logic_error
        {
        public:
            MatrixInversionException(const std::string& __arg) :
            std::logic_error(__arg) {
            }

        };
        
            template<typename NT, unsigned long N, unsigned long M>
            class MathMatrix;
            
            template<typename NT, unsigned long N, unsigned long M>
            class MatrixHelper;
            
            template<typename NT, unsigned long N, unsigned long M>
            class MathMatrix
            {
                public:
                    typedef MathVector<NT, M, true> MatrixColumn;
                private:
                    MatrixColumn columns_[N];
                protected:
                    template<typename NT2, unsigned long N2, unsigned long M2>
                    friend class MathMatrix;
                    
                    class MatrixRow
                    {
                        MathMatrix<NT,N,M>& parent;
                        unsigned long y;
                    public:
                        virtual ~MatrixRow() {
                        }
                        
                        MatrixRow(MathMatrix<NT, N, M>& parent, unsigned long y) :
                            parent(parent), y(y) {
                        }
                        
                        NT& operator[](unsigned long index) {
                            return parent[index][y];
                        }

                        
                        MatrixRow& operator*=(NT right) {
                            for(unsigned long x = 0; x < N; x++)
                            {
                                parent[x][y] *= right;
                            }
                            return *this;
                        }

                        template<typename IV>
                        MatrixRow& operator+=(const IV& right) {
                            for(unsigned long x = 0; x < N; x++)
                            {
                                parent[x][y] += right[x];
                            }
                            return *this;
                        }

                        template<typename IV>
                        MatrixRow& operator-=(const IV& right) {
                            for(unsigned long x = 0; x < N; x++)
                            {
                                parent[x][y] -= right[x];
                            }
                            return *this;
                        }

                        MatrixRow& operator/=(NT right) {
                            for(unsigned long x = 0; x < N; x++)
                            {
                                parent[x][y] /= right;
                            }
                            return *this;
                        }
                        
                        MathVector<NT,N,false> asVector()
                        {
                            MathVector<NT,N,false> ret;
                            for(unsigned long x = 0; x < N; x++)
                                ret[x] = (*this)[x];
                            return ret;
                        }
                    };
                    
                    MatrixRow row(unsigned long y)
                    {
                        return MatrixRow(*this,y);
                    }
                    
                    void swapRows(unsigned long y1, unsigned long y2)
                    {
                        for(unsigned long x = 0; x < N; x++)
                        {
                            std::swap(columns_[x][y1], columns_[x][y2]);
                        }
                    }
                    
                public:
                    MathMatrix(const MatrixColumn& column)
                    {
                        for(auto i = 0; i < N; i++)
                        {
                            columns_[i] = column;
                        }
                    }
                    
                    MathMatrix(const MathVector<NT, N, false>& row)
                    {
                        for(auto y = 0; y < M; y++)
                        {
                            for(auto x = 0; x < N; x++)
                            {
                                columns_[x][y] = row[x];
                            }
                        }
                    }
                    
                    MathMatrix(const std::initializer_list<NT>& values)
                    {
                        unsigned long minInit = std::min(N*M,values.size()); 
                        unsigned long i = 0;
                        
                        for(auto it = values.begin(); it != values.end(); it++)
                        {
                            columns_[i % N][i / N] = *it;
                            i++;
                            if(i >= minInit)
                                break;
                        }
                    }
                    
                    MathMatrix(NT init) 
                    {
                        for(int i = 0; i < N; i++)
                            columns_[i] = init;
                    }

                    
                    MathMatrix()
                    {
                    }
                    
                    static MathMatrix<NT, N, M> identity()
                    {
                        MathMatrix<NT, N, M> ret;
                        for(unsigned long x = 0; x < N; x++)
                        {
                            for(unsigned long y = 0; y < M; y++)
                            {
                                ret[x][y] = x==y ? 1 : 0;
                            }
                        }
                        return ret;
                    }
                    
                    MatrixColumn& operator[](unsigned long columnNumber) {
                        return columns_[columnNumber];
                    }
                    
                    const MatrixColumn& operator[](unsigned long index) const {
                        return const_cast<MathMatrix<NT, N, M>&> (*this)[index];
                    }

                    
                    
                    MathMatrix<NT, N, M>& operator*=(NT right) {
                        for(auto& row : columns_)
                        {
                            row *= right;
                        }
                        return *this;
                    }
                    
                    MathMatrix<NT, N, M>& operator/=(NT right) {
                        for(auto& row : columns_)
                        {
                            row /= right;
                        }
                        return *this;
                    }


                    MathMatrix<NT, N, M>& operator+=(const MathMatrix<NT, N, M>& right) {
                        for(unsigned long i = 0; i < N; i++)
                            columns_[i] += right.columns_[i];
                        return *this;
                    }

                    MathMatrix<NT, N, M>& operator-=(const MathMatrix<NT, N, M>& right) {
                        for(unsigned long i = 0; i < N; i++)
                            columns_[i] -= right.columns_[i];
                        return *this;
                    }
                    
                    MathMatrix<NT, M, N> transposition() const
                    {
                        MathMatrix<NT, M, N> ret;
                        for(long unsigned x = 0; x < N; x++)
                            for(long unsigned y = 0; y < M; y++)
                            {
                                ret[y][x] = (*this)[x][y];
                            }
                        return ret;
                    }
                    
                    MathMatrix<NT, M, N> diag() const
                    {
                        MathMatrix<NT, M, N> ret;
                        for(long unsigned x = 0; x < N; x++)
                            for(long unsigned y = 0; y < M; y++)
                            {
                                if(x == y)
                                    ret[x][y] = (*this)[x][y];
                                else
                                    ret[x][y] = 0;
                            }
                        return ret;
                    }
                    
                    template<unsigned long N2> 
                    MathMatrix<NT, N+N2, M> concat(const MathMatrix<NT, N2, M>& mat) const
                    {
                        MathMatrix<NT,N+N2, M> ret;
                        for(unsigned long x = 0; x < N; x++)
                        {
                            for(unsigned long y = 0; y < M; y++)
                            {
                                ret[x][y] = (*this)[x][y];
                            }
                        }
                        
                        for(unsigned long x = 0; x < N2; x++)
                        {
                            for(unsigned long y = 0; y < M; y++)
                            {
                                ret[N+x][y] = mat[x][y];
                            }
                        }
                        return ret;
                    }
                    
                    template<unsigned long SX, unsigned long SY, unsigned long LX = N - SX, unsigned long LY = M - SY>
                    MathMatrix<NT, LX, LY> submatrix() const
                    {
                        MathMatrix<NT, LX, LY> ret;
                        for(unsigned long x = 0; x < LX; x++)
                        {
                            for(unsigned long y = 0; y < LY; y++)
                            {
                                ret[x][y] = (*this)[SX+x][SY+y];
                            }
                        }
                        return ret;
                    }
                    
                    MathMatrix<NT,N - 1,M - 1> matrixMinor(unsigned long X, unsigned long Y) const
                    {
                        MathMatrix<NT,N - 1,M - 1> ret;
                        
                        unsigned long sx = 0, sy = 0, dx = 0, dy = 0;
                        for(sx = 0; sx < N; sx++ )
                        {
                            if(sx != X)
                            {
                                dy = 0;
                                for(sy = 0; sy < M; sy++ )
                                {
                                    if(sy != Y)
                                    {
                                        ret[dx][dy] = (*this)[sx][sy];
                                        dy++;
                                    }
                                }
                                dx++;
                            }
                        }
                        return ret;
                    }
                    
                    NT det() const
                    {
                        return MatrixHelper<NT,N,M>::MatrixDeterminant(*this);
                    }
                    
                    MathMatrix<NT, N, M> inverse() const
                    {
                        if(N != M)
                            throw std::logic_error("It is possible to inverse only square marix.");
                        

                        if(det() == 0)
                            throw MatrixInversionException("Matrix is singular - cannot be reversed");

                        
                        auto concatenated = this->concat(MathMatrix<NT, N, M>::identity());
                        
                        for(unsigned long y = 0; y < M; y++)
                        {
                            do
                            {
                                unsigned long startX;
                                for(unsigned long x = y; x < N; x++)
                                    if(concatenated[x][y] != 0)
                                    {
                                        startX = x;
                                        break;
                                    }
                                
                                if(startX != y)
                                {
                                    concatenated.swapRows(y,startX);
                                }
                            }
                            while(concatenated[y][y] == 0);
                        }
                        
                        for(unsigned long y = 0; y < M; y++)
                        {
                            unsigned long startX = y;
                            
                            NT v = concatenated[startX][y];
                            concatenated.row(y) /= v;
                            
                            for(unsigned long y2 = 0; y2 < M; y2++)
                            {
                                if(y == y2)
                                    continue;
                                auto extractedRow = concatenated.row(y).asVector();
                                extractedRow *= concatenated[startX][y2];
                                concatenated.row(y2) -= extractedRow;
                            }
                        }
                        
                        
                        return concatenated.template submatrix<N,0>();
                    }
                    
                    template<class Archive>
                    void serialize(Archive & ar, const unsigned int version)
                    {
                        ar & columns_;
                    }
                    
                    static MatrixColumn gaussJordanSolve(const MathMatrix<NT, N, M>& A, const MatrixColumn& B)
                    {
                        if(N != M)
                            throw std::logic_error("Gauss - Jordan method is implmemented only for square matrixes");
                        

                        if(A.det() == 0)
                            throw MatrixInversionException("A matrix is singular - solution cannod be found");
                        
                        auto concatenated = A.concat(MathMatrix<NT,1,M>(B));
                        
                        for(unsigned long y = 0; y < M; y++)
                        {
                            do
                            {
                                unsigned long startX;
                                for(unsigned long x = y; x < N; x++)
                                    if(concatenated[x][y] != 0)
                                    {
                                        startX = x;
                                        break;
                                    }
                                
                                if(startX != y)
                                {
                                    concatenated.swapRows(y,startX);
                                }
                            }
                            while(concatenated[y][y] == 0);
                        }
                        
                        for(unsigned long y = 0; y < M; y++)
                        {
                            unsigned long startX = y;
                            
                            NT v = concatenated[startX][y];
                            concatenated.row(y) /= v;
                            
                            for(unsigned long y2 = 0; y2 < M; y2++)
                            {
                                if(y == y2)
                                    continue;
                                auto extractedRow = concatenated.row(y).asVector();
                                extractedRow *= concatenated[startX][y2];
                                concatenated.row(y2) -= extractedRow;
                            }
                        }
                        
                        return concatenated.template submatrix<N,0>().columns_[0];
                    }
            };
            
            template<typename NT, unsigned long N, unsigned long M>
            MathMatrix<NT, N, M> operator*(const MathMatrix<NT, N, M>& left, NT right) {
                MathMatrix<NT, N, M> ret = left;
                ret *= right;
                return ret;
            }
            
            template<typename NT, unsigned long N, unsigned long M>
            MathMatrix<NT, N, M> operator/(const MathMatrix<NT, N, M>& left, NT right) {
                MathMatrix<NT, N, M> ret = left;
                ret /= right;
                return ret;
            }
            
            template<typename NT, unsigned long N, unsigned long M>
            MathMatrix<NT, N, M> operator+(const MathMatrix<NT, N, M>& left, const MathMatrix<NT, N, M>& right) {
                MathMatrix<NT, N, M> ret = left;
                ret += right;
                return ret;
            }
            
            template<typename NT, unsigned long N, unsigned long M>
            MathMatrix<NT, N, M> operator-(const MathMatrix<NT, N, M>& left, const MathMatrix<NT, N, M>& right) {
                MathMatrix<NT, N, M> ret = left;
                ret -= right;
                return ret;
            }
            
            template<typename NT, unsigned long N, unsigned long M, unsigned long K>
            MathMatrix<NT, K, M> operator*(const MathMatrix<NT, N, M>& left, const MathMatrix<NT, K, N>& right) {
                MathMatrix<NT, K, M> ret;
                for(unsigned long y = 0; y < M; y++)
                {
                    for(unsigned long x = 0; x < K; x++)
                    {
                        NT sum = 0;
                        for(unsigned long n = 0; n < N; n++)
                        {
                            sum += left[n][y] * right[x][n];
                        }
                        ret[x][y] = sum;
                    }
                }
                
                return ret;
            }
            
            template<typename NT, unsigned long N, unsigned long M>
            MathVector<NT, N> operator*(const MathMatrix<NT, N, M>& left, const MathVector<NT, N>& right) {
                MathVector<NT, N> ret;
                for(unsigned long y = 0; y < N; y++)
                {
                    NT sum = 0;
                    for(unsigned long n = 0; n < N; n++)
                    {
                        sum += left[n][y] * right[n];
                    }
                    ret[y] = sum;
                }
                
                return ret;
            }
            
            template<typename NT, unsigned long N, unsigned long M>
            MathVector<NT, N, false> operator*(const MathVector<NT, M, false>& left, const MathMatrix<NT, N, M>& right) {
                MathVector<NT, N, false> ret;
                for(unsigned long x = 0; x < M; x++)
                {
                    NT sum = 0;
                    for(unsigned long n = 0; n < M; n++)
                    {
                        sum += left[n] * right[x][n];
                    }
                    ret[x] = sum;
                }
                
                return ret;
            }
            
            template<typename NT, unsigned long N>
            NT operator*(const MathVector<NT, N, false>& left, const MathVector<NT, N, true>& right) 
            {
                NT sum = 0;
                for(unsigned long i = 0; i < N; i++)
                    sum += left[i] * right[i];
                return sum;
            }
            
            template<typename NT, unsigned long N>
            MathMatrix<NT,N,N> operator*(const MathVector<NT, N, true>& left, const MathVector<NT, N, false>& right) 
            {
                MathMatrix<NT,N,N> ret;
                for(unsigned long x=0; x<N; x++)
                    for(unsigned long y=0; y<N; y++)
                        ret[x][y] = left[x] * right[y];
                return ret;
            }
            
            
            template<typename NT, unsigned long N, unsigned long M>
            class MatrixHelper
            {
                public:
                static NT MatrixDeterminant(const MathMatrix<NT,N,M>& mat)
                {  
                    int sign = 1;
                    NT sum = 0;
                    for(long unsigned l = 0; l < N; l++)
                    {
                        auto minorMat = mat.matrixMinor(0,l);
                        auto minorMatDet = minorMat.det();
                        sum += sign * mat[0][l]*minorMatDet;
                        sign = -sign;
                    }

                    return sum;
                }
            };
            
            template<typename NT>
            class MatrixHelper<NT,1,1>
            {
                public:
                
                static NT MatrixDeterminant(const MathMatrix<NT,1,1>& mat)
                {
                    return mat[0][0];
                }
            };
    }
    
}
}

#endif	/* MATHMATRIX_H */

