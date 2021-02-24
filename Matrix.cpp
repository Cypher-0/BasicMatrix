#include "Matrix.hpp"

#include <iostream>

#include <algorithm>

namespace mat
{

namespace test
{

void matrix()
{
    mat::Matrix<int> mat{{1,5,-3},{3,4,-1},{6,8,-2}};
    mat::Matrix<int> mat1{mat.m(),mat.n(),1};

    printSep("Transpose");//---------------------------------------------

    dispMatrix("Base matrix",mat);
    dispMatrix("Transposed",t(mat));

    printSep("Addition");//---------------------------------------------

    dispMatrix("Base matrix M0",mat);
    dispMatrix("Base matrix M1",mat1);
    dispMatrix("M0 + M1",mat+mat1);

    printSep("Multiplication");//---------------------------------------------

    dispMatrix("Base matrix M0",mat);
    dispMatrix("Base matrix M1",mat1);
    dispMatrix("M0 * M1",mat*mat1);
    printSep();

    dispMatrix("Base matrix M0",mat);
    dispMatrix("M0 * 2",mat*2);
    printSep();

    dispMatrix("Base matrix M0",mat);
    mat *= 2;
    dispMatrix("M0 *= 2",mat);


    printSep("Division");//---------------------------------------------

    dispMatrix("Base matrix M0",mat);
    dispMatrix("M0 / 2",mat/2);
    printSep();

    dispMatrix("Base matrix M0",mat);
    mat /= 2;
    dispMatrix("M0 /= 2",mat);


    printSep("Identity");//---------------------------------------------

    dispMatrix("Id_2",identity<int>(2));
    dispMatrix("Id_3",identity<int>(3));
    dispMatrix("Id_4",identity<int>(4));
    printSep();

    dispMatrix("Base matrix M0",mat);
    dispMatrix("M0 * Id_3",mat*identity<int>(3));


    printSep("Generate random");//---------------------------------------------

    dispMatrix("Normal randomized",generateNormalRandom(3,3));


    printSep("Split y");//---------------------------------------------

    dispMatrix("Base matrix M0",mat);
    auto [part0,part1]{mat.split(1)};
    dispMatrix("Part0 (split at 1)",part0);
    dispMatrix("Part1 (split at 1)",part1);




    printSep(" Temporary ");//--------------------------------------------- TEMP TESTs
    dispMatrix("Base",mat1);
}

}//namespace test

}//namespace mat

