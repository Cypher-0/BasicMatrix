#ifndef MAT_MATRIX_HPP
#define MAT_MATRIX_HPP


#include <vector>
#include <string>
#include <tuple>
#include <optional>

#include <random>
#include <algorithm>

#include <cassert>
#include <iostream>

namespace mat
{

template<typename T = double>
class Matrix
{
public:
    using size_type = typename std::vector<T>::size_type;
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;
    //----------------------------------------------------- SPECIALs

    Matrix(size_type m,size_type n,const T& defVal=T{}) : m_datas{},m_m{m},m_n{n}
    {
        m_datas.resize(m_m*m_n,defVal);
    }

    Matrix(std::initializer_list<std::vector<T>> datas) : m_datas{},m_m{},m_n{}
    {
        std::vector<std::vector<T>> temp{std::move(datas)};
        size_type baseSize{};
        if(size(temp) > 0)
            baseSize = size(temp[0]);

        m_m = size(temp);
        m_n = baseSize;

        size_type index{};
        m_datas.resize(m_m*m_n);

        for(const auto& e : temp)
        {
            assert(size(e) == baseSize && "Every row of the matrix has to have the same size");

            for(const auto& el : e)
            {
                m_datas[index] = std::move(el);
                ++index;
            }
        }
    }

    /*
    Matrix(std::vector<std::vector<T>> datas) : m_datas{std::move(datas)},m_m{},m_n{}
    {
        using std::to_string;
        using std::size;
        size_type baseSize{};
        if(size(m_datas) > 0)
            baseSize = size(m_datas[0]);

        for(const auto& e : m_datas)
        {
            assert(size(e) == baseSize && "Every row of the matrix has to have the same size");
        }
        m_m = size(m_datas);
        m_n = baseSize;
    }*/

    Matrix() : Matrix(0,0)
    {
        //default constructor
    }

    ~Matrix()
    {

    }

////////////////////////////////////
////
////    Accessors
////
////////////////////////////////////

const std::vector<T>& datas() const {
    return m_datas;
}

size_type m() const{ //return rows number
    return m_m;
}

size_type n() const{ //return column number
    return m_n;
}

const T& at(size_type i,size_type j) const{ //return const ref to element at
    return m_datas[j*m_n+i];
}

T& at(size_type i,size_type j){ //return reference to element at
    //assert(i < m_n && j < m_m);
    return m_datas[j*m_n+i];
}

    //raw array access

const T& at(size_type i) const{ //return const ref to element at
    return m_datas[i];
}

T& at(size_type i){ //return reference to element at
    return m_datas[i];
}

////////////////////////////////////
////
////    Tools (split)
////
////////////////////////////////////

std::tuple<Matrix<T>,Matrix<T>> split(typename Matrix<T>::size_type index, bool x = false) //x represents the axis on which we will split (x = false imply split on y)
{
    return (x)?M_splitx(index):M_splity(index);
}


////////////////////////////////////
////
////    Operators
////
////////////////////////////////////

    //Additions . substractions

Matrix<T>& operator+=(const Matrix<T>& mat){
    *this = *this+mat;
    return *this;
}

Matrix<T>& operator-=(const Matrix<T>& mat){
    *this = *this-mat;
    return *this;
}

    //Multiplications . divisions

Matrix<T>& operator*=(const Matrix<T>& mat){
    *this = *this * mat;
    return *this;
}

Matrix<T>& operator*=(const T& n){
    *this = *this * n;
    return *this;
}

Matrix<T>& operator/=(const T& n){
    *this = *this / n;
    return *this;
}

////////////////////////////////////
////
////    Iterators
////
////////////////////////////////////

iterator begin(){
    return std::begin(m_datas);
}

const_iterator cbegin() const {
    return std::cbegin(m_datas);
}

iterator end(){
    return std::end(m_datas);
}

const_iterator cend() const {
    return std::cend(m_datas);
}

//********************
//  Private Members
//********************

private:
    std::vector<T> m_datas;

    size_type m_m;//rows count
    size_type m_n;//columns count

    ////////////////////////////////////
    ////
    ////    Tools (split)
    ////
    ////////////////////////////////////

std::tuple<Matrix<T>,Matrix<T>> M_splity(typename Matrix<T>::size_type rowIndex)
{
    Matrix<T> out0{rowIndex,n()};
    Matrix<T> out1{m()-rowIndex,n()};

    auto index{rowIndex*n()};

    std::copy(this->cbegin(),this->cbegin()+index,std::begin(out0.m_datas));
    std::copy(this->cbegin()+index,this->cend(),std::begin(out1.m_datas));

    return {out0,out1};
}

std::tuple<Matrix<T>,Matrix<T>> M_splitx(typename Matrix<T>::size_type colIndex)
{
    assert("" == "WORK IN PROGRESS");
    return {{},{}};
}

};//----------------------------------------------------------------------------------------------------------------------------------------

////################################
////################################
////
////    operators
////
////################################
////################################


    //----------------------------
    //
    //  addition . substraction
    //
    //----------------------------

template<typename T>
Matrix<T> operator+(const Matrix<T>& m0,const Matrix<T>& m1)
{
    assert(m0.m() == m1.m() && m0.n() == m1.n());

    Matrix<T> out{m0.m(),m0.n()};
    for(auto i{m0.m()};i > 0; --i)
    {
        for(auto j{m0.n()}; j > 0; --j)
        {
            auto y{i-1};
            auto x{j-1};
            out.at(x,y) = m0.at(x,y)+m1.at(x,y);
        }
    }
    return out;
}

template<typename T>
Matrix<T> operator-(const Matrix<T>& m0,const Matrix<T>& m1)
{
    assert(m0.m() == m1.m() && m0.n() == m1.n());

    Matrix<T> out{m0.m(),m0.n()};
    for(auto i{m0.m()};i > 0; --i)
    {
        for(auto j{m0.n()}; j > 0; --j)
        {
            auto y{i-1};
            auto x{j-1};
            out.at(x,y) = m0.at(x,y)-m1.at(x,y);
        }
    }
    return out;
}

    //----------------------------
    //
    //      multiplication
    //
    //----------------------------

template<typename T>
Matrix<T> operator*(const Matrix<T>& a,const Matrix<T>& b)
{
    assert(b.m() == a.n());

    Matrix<T> out{a.m(),b.n()};//init an empty new matrix with the appropriates dimensions
    using size_type = typename Matrix<T>::size_type;
    for(size_type i{}; i < a.m(); ++i)
    {
        for(size_type j{}; j < b.n(); ++j)
        {
            for(size_type k{}; k < a.n(); ++k)
            {
                out.at(j,i) += a.at(k,i)*b.at(j,k);
            }
        }
    }

    return out;
}

template<typename T>
Matrix<T> multElems(const Matrix<T>& a,const Matrix<T>& b)//out[i,j] = a[i,j]*b[i,j]
{
    assert(a.m() == b.m() && a.n() == b.n());

    Matrix<T> out{a.m(),b.n()};//init an empty new matrix with the appropriates dimensions
    using size_type = typename Matrix<T>::size_type;
    for(size_type i{}; i < a.m(); ++i)
    {
        out.at(i) = a.at(i)*b.at(i);
    }

    return out;
}

    //----------------------------
    //
    //      Scalars operators
    //
    //----------------------------

template<typename T>
Matrix<T> operator*(const Matrix<T>& m0,const T& n)
{
    Matrix<T> out{m0.m(),m0.n()};

    for(auto i{out.m()};i > 0; --i)
    {
        for(auto j{out.n()}; j > 0; --j)
        {
            auto y{i-1};
            auto x{j-1};

            out.at(x,y) = m0.at(x,y)*n;
        }
    }
    return out;
}
template<typename T>
Matrix<T> operator*(const T& n,const Matrix<T>& m0){
    return m0*n;
}

template<typename T>
Matrix<T> operator/(const Matrix<T>& m0,const T& n)
{
    Matrix<T> out{m0.m(),m0.n()};

    for(auto i{out.m()};i > 0; --i)
    {
        for(auto j{out.n()}; j > 0; --j)
        {
            auto y{i-1};
            auto x{j-1};
            auto temp{m0.at(x,y)/n};
            out.at(x,y) = std::move(temp);
        }
    }
    return out;
}

////################################
////################################
////
////    Maths
////
////################################
////################################

    //-----------------------
    //
    //  TRANSPOSE
    //
    //-----------------------

template<typename T>
Matrix<T> swapLinesAndRows(const Matrix<T>& mat) //make a 3x2 matrix become a 2x3 matrix == transpose on a square matrix
{
    Matrix<T> out{mat.n(),mat.m()};

    for(typename Matrix<T>::size_type i{}; i < mat.n();++i)
    {
        for(typename Matrix<T>::size_type j{}; j < mat.m();++j)
        {
            out.at(j,i) = mat.at(i,j);
        }
    }
    return out;
}

template<typename T>
Matrix<T> transpose(const Matrix<T>& mat)
{
    assert (mat.m() == mat.n());
    return swapLinesAndRows(mat);
}

template<typename T>
Matrix<T> t(const Matrix<T> &mat)
{
    return transpose(mat);
}

    //-----------------------
    //
    //  Identity
    //
    //-----------------------

template<typename T>
Matrix<T> identity(typename Matrix<T>::size_type dim)
{
    Matrix<T> out{dim,dim};
    for(typename Matrix<T>::size_type i{}; i < dim;++i)
    {
        out.at(i,i) = static_cast<T>(1);
    }
    return out;
}

    //-----------------------
    //
    //  Maxs
    //
    //-----------------------

template<typename T>
Matrix<T> max_columns(const Matrix<T>& mat)//search the max value in each columns and return them as an array
{
    if(mat.m() < 1)
        return Matrix<T>{};

    Matrix<T> out{1,mat.n()};

    for(typename Matrix<T>::size_type i{}; i < mat.n(); ++i)//go through all columns of the matrix
    {
        T maxVal{mat.at(i,0)};

        for(typename Matrix<T>::size_type j{1}; j < mat.m(); ++j)
        {
            maxVal = std::max(maxVal,mat.at(i,j));
        }

        out.at(i,0) = maxVal;
    }

    return out;
}

    //-----------------------
    //
    //  Normalize
    //
    //-----------------------

template<typename T>
Matrix<T> normalizeByColumns(const Matrix<T>& mat)//make sure each column of the matrix is normalized
                                                  //each value of each column get divided by maximum value of the column
{
    auto maxVals{max_columns(mat)};
    Matrix<T> out{mat.m(),mat.n()};

    for(typename Matrix<T>::size_type i{}; i < mat.n();++i)
    {
        auto maxVal{maxVals.at(i,0)};
        for(typename Matrix<T>::size_type j{}; j < mat.m();++j)
        {
            out.at(i,j) = mat.at(i,j)/maxVal;
        }
    }

    return out;
}


////################################
////################################
////
////    Iterators
////
////################################
////################################

template<typename T>
typename Matrix<T>::iterator begin(Matrix<T>& mat){
    return mat.begin();
}

template<typename T>
typename Matrix<T>::const_iterator cbegin(const Matrix<T>& mat){
    return mat.cbegin();
}


template<typename T>
typename Matrix<T>::iterator end(Matrix<T>& mat){
    return mat.end();
}

template<typename T>
typename Matrix<T>::const_iterator cend(const Matrix<T>& mat){
    return mat.cend();
}

////################################
////################################
////
////    Random generation
////
////################################
////################################


inline
Matrix<double> generateRandom(typename Matrix<double>::size_type m,typename Matrix<double>::size_type n,double min, double max)
{
    static std::random_device rd;  //Will be used to obtain a seed for the random number engine
    static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> dis(min, max);

    Matrix<double> mat{m,n};

    auto genRand{
    [&](auto& e)
    {
        e = dis(gen);
    }
    };

    std::for_each(begin(mat),end(mat),genRand);
    return mat;
}


inline
Matrix<double> generateNormalRandom(typename Matrix<double>::size_type m,typename Matrix<double>::size_type n,double min=0.0, double max=1.0)
{
    static std::random_device rd;  //Will be used to obtain a seed for the random number engine
    static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd() mt19937
    std::normal_distribution<double> dis(min, max);

    Matrix<double> mat{m,n};

    auto genRand{
    [&](auto& e)
    {
        e = dis(gen);
    }
    };
    std::for_each(begin(mat),end(mat),genRand);
    return mat;
}

inline
Matrix<int> generateRandom(typename Matrix<int>::size_type m,typename Matrix<int>::size_type n,int min, int max)
{
    static std::random_device rd;  //Will be used to obtain a seed for the random number engine
    static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<int> dis(min, max);

    Matrix<int> mat{m,n};

    auto genRand{
    [&](auto& e)
    {
        e = dis(gen);
    }
    };
    std::for_each(begin(mat),end(mat),genRand);
    return mat;
}


////################################
////################################
////
////    String and formatting
////
////################################
////################################

template<typename T>
std::string to_string(const Matrix<T>& mat)
{
    using std::to_string;

    using size_type = typename Matrix<T>::size_type;

    std::string out{"["};

    for(size_type j{}; j < mat.m();++j)
    {
        out += "[";
        for(size_type i{}; i < mat.n();++i)
        {
            out += to_string(mat.at(i,j))+",";
        }
        if(mat.n() > 0)
        {
            out.pop_back();
        }
        out += "]\n";
    }

    if(size(mat.datas())>0)//if there was at least 1 element
        out.pop_back(); // remove the last char (`\n`)
    out += "]";
    return out;
}

template<typename T>
std::ostream& operator<<(std::ostream& stream,const Matrix<T>& mat)
{
    return stream << to_string(mat);
}


////===============================================================================================
////
////                                TESTs namespace
////
////===============================================================================================

namespace test
{
template<typename T>
void dispMatrix(const std::string& name,const Matrix<T>& mat)
{
    std::cout << name <<" : \n";
    std::cout << to_string(mat)<<"\n\n";
}

inline
void printSep(const std::string& title="")
{
    std::cout << "\t-----------------------------"+title+"----------------------------\n";
}

void matrix();

}//namespace test

}//namespace mat


#endif //MAT_MATRIX_HPP
