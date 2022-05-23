/* 
 * File:   MathVector.h
 * Author: tombieli
 *
 * Created on 6 November 2015, 15:19
 */
#ifndef MATHVECTOR_H
#define	MATHVECTOR_H


#include <cmath>
#include <cstring>
#include <iterator>
#include <stdexcept>
#include <functional>

namespace ksg { 

namespace utils {

namespace math
{
            template<typename NT, unsigned long S, bool V = true>
            class MathVector
            {
                friend class MathVector<NT,S,!V>;
                NT content_[S];
            public:
                typedef NT number_type;
                typedef NT* iterator;
                typedef const NT* const_iterator;
                typedef std::reverse_iterator<iterator> reverse_iterator;
                typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
                static const unsigned long size = S;
                static const bool isVertical = V;
                MathVector(NT init)
                {
                    for(unsigned long i=0; i < S; i++)
                    {
                        content_[i] = init; 
                    }
                }
                
                template<bool OV>
                MathVector(const MathVector<NT, S, OV>& right) {
                    std::copy(right.content_, right.content_+S, content_);
                }
                
                MathVector(std::initializer_list<NT> list) {
                    if(list.size() != S)
                        throw std::runtime_error("Bad initializer list size");
                    
                    unsigned long i = 0;
                    for(auto e : list)
                    {
                        content_[i] = e;
                        i++;
                    }
                }
                
                MathVector()
                    : MathVector(0)
                {
                }
                
                iterator begin()
                {
                    return content_;
                }
                
                const_iterator begin() const
                {
                    return content_;
                }
                
                iterator end()
                {
                    return content_ + size;
                }
                
                const_iterator end() const
                {
                    return content_ + size;
                }
                
                reverse_iterator rbegin()
                {
                    return reverse_iterator(&content_[size - 1]);
                }
                
                const_reverse_iterator rbegin() const
                {
                    return const_reverse_iterator(&content_[size - 1]);
                }
                
                reverse_iterator rend()
                {
                    return reverse_iterator(content_ - 1);
                }
                
                const_reverse_iterator rend() const
                {
                    return const_reverse_iterator(content_ - 1);
                }
                
                NT& operator[](unsigned long index) {
                    return content_[index];
                }

                const NT& operator[](unsigned long index) const {
                    return const_cast<MathVector<NT,S,V>&> (*this)[index];
                }                
                
                MathVector<NT,S,V>& operator=(const MathVector<NT,S,V>& right) {
                    if (this == &right)
                        return *this; 
                    
                    std::copy(right.content_, right.content_ + S, this->content_);
                    
                    return *this;
                }
                
                NT abs()
                {
                    NT sum = 0;
                    for(int i = 0; i < S; i++)
                    {
                        auto v = (*this)[i];
                        sum += v*v;
                    }
                    return sqrt(sum);
                }

                template<typename NT2>
                MathVector<NT,S,V>& operator*=(NT2 right) {
                    for(int i=0; i < S; i++)
                    {
                        (*this)[i] *= right;
                    }
                    return *this;
                }

                MathVector<NT,S,V>& operator+=(const MathVector<NT,S,V>& right) {
                    for(int i=0; i < S; i++)
                    {
                        (*this)[i] += right[i];
                    }
                    return *this;
                }

                MathVector<NT,S,V>& operator-=(const MathVector<NT,S,V>& right) {
                    for(int i=0; i < S; i++)
                    {
                        (*this)[i] -= right[i];
                    }
                    return *this;
                }


                template<typename NT2>
                MathVector<NT,S,V>& operator/=(const NT2& right) {
                    for(int i=0; i < S; i++)
                    {
                        (*this)[i] /= right;
                    }
                    return *this;
                }

                MathVector<NT,S,V> operator-()
                {
                    MathVector<NT,S,V> ret((NT)0.0);
                    ret -= *this;
                    return ret;
                }
                
                MathVector<NT,S,V>& operator&=(const MathVector<NT,S,V>& right) {
                    for(unsigned long i =0; i < S; i++)
                        (*this)[i] *= right[i];
                    return *this;
                }
                
                MathVector<NT,S,V>& genericInPlaceOperation(std::function<NT(NT)> unary_operator) {
                    for(unsigned long i =0; i < S; i++)
                        (*this)[i] = unary_operator((*this)[i]);
                    return *this;
                }
                
                MathVector<NT,S,V>& genericInPlaceOperation(const MathVector<NT,S,V>& right, std::function<NT(NT,NT)> binary_operator) {
                    for(unsigned long i =0; i < S; i++)
                        (*this)[i] = binary_operator((*this)[i], right);
                    return *this;
                }
                

                
                NT dot(const MathVector<NT,S,V>& right)
                {
                    NT sum = 0;
                    for(unsigned long i =0; i < S; i++)
                        sum += (*this)[i] * right[i];
                       
                    return sum;
                }
                
                MathVector<NT,S,V> clamp(const MathVector<NT,S,V>& min, const MathVector<NT,S,V>& max)
                {
                    MathVector<NT,S,V> ret;
                    
                    for(unsigned long i=0; i < S; i++)
                    {
                        NT value = (*this)[i];
                        if(value < min[i])
                            ret[i] = min[i];
                        else if(value > max[i])
                            ret[i] = max[i];
                        else 
                            ret[i] = value;
                    }
                    
                    return ret;
                }
                
                MathVector<NT,S,V> clampRescaleABS(const MathVector<NT,S,V>& bounds)
                {
                    MathVector<NT,S,V> ret(*this);
                    
                    for(unsigned long i=0; i < S; i++)
                    {
                        NT value = ret[i];
                        NT absValue = value < 0 ? -value : value;
                        if(absValue > bounds[i])
                        {
                            NT coef = bounds[i] / absValue;
                            ret *= coef;
                        }
                    }
                    
                    return ret;
                }
                
                MathVector<NT,S,!V> transposition()
                {
                    return MathVector<NT,S,!V>(*this);
                }
                
                MathVector<NT,S-1,V> vectorMinor(unsigned long index) const 
                {
                    MathVector<NT,S-1,V> ret;
                    for(unsigned long i=0; i < index; i++)
                        ret[i] = (*this)[i];
                    for(unsigned long i=index + 1; i < S; i++)
                        ret[i-1] = (*this)[i];
                        
                    return ret;
                }
                
                MathVector<NT,S+1,V> vectorInserValue(unsigned long index, NT value) const
                {
                    MathVector<NT,S+1,V> ret;
                    for(unsigned long i=0; i < index; i++)
                        ret[i] = (*this)[i];
                    
                        ret[index] = value;
                        
                    for(unsigned long i=index; i < S; i++)
                        ret[i+1] = (*this)[i];
                        
                    return ret;
                }
                
                template<class Archive>
                void serialize(Archive & ar, const unsigned int version)
                {
                    ar & content_;
                }
                

            };
            
            template<typename NT, unsigned long S, bool V>
            MathVector<NT,S,V> genericOperation(const MathVector<NT,S,V>& left, const MathVector<NT,S,V>& right, std::function<NT(NT,NT)> binary_operator) {
                MathVector<NT,S,V> result(left); // Make a copy of myself.
                result.genericInPlaceOperation(right, binary_operator); // Reuse compound assignment
                return result;
            }
            
            template<typename NT, unsigned long S,bool V, typename NT2>
            MathVector<NT,S,V> operator*(const MathVector<NT,S,V>& left, NT2 right) {
                MathVector<NT,S,V> result(left); // Make a copy of myself.
                result *= right; // Reuse compound assignment
                return result;
            }

            template<typename NT, unsigned long S,bool V>
            MathVector<NT,S,V> operator+(const MathVector<NT,S,V>& left,const MathVector<NT,S,V>& right) {
                MathVector<NT,S,V> result(left); // Make a copy of myself.
                result += right; // Reuse compound assignment
                return result;
            }

            template<typename NT, unsigned long S,bool V>
            MathVector<NT,S,V> operator-(const MathVector<NT,S,V>& left,const MathVector<NT,S,V>& right) {
                MathVector<NT,S,V> result(left); // Make a copy of myself.
                result -= right; // Reuse compound assignment
                return result;
            }

            template<typename NT, unsigned long S,bool V, typename NT2>
            MathVector<NT,S,V> operator/(const MathVector<NT,S,V>& left, NT2 right) {
                MathVector<NT,S,V> result(left); // Make a copy of myself.
                result /= right; // Reuse compound assignment
                return result;
            }
            
            template<typename NT, unsigned long S,bool V>
            MathVector<NT,S,V> operator&(const MathVector<NT,S,V>& left, const MathVector<NT,S,V>& right){
                MathVector<NT,S,V> result(left); // Make a copy of myself.
                result &= right; // Reuse compound assignment
                return result;
            }
            
            
            
}
    
}
}

#endif	/* MATHVECTOR_H */

