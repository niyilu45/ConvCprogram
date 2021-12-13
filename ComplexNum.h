#ifndef COMPLEXNUM_H
#define COMPLEXNUM_H

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef COMPLEXNUM
#define COMPLEXNUM
typedef struct ComplexNumStru
{
    double re;
    double im;
} ComplexNumStru;
typedef struct ComplexNumStru ComplexNum;
#endif

ComplexNum ComplexNumReset(void);
ComplexNum ComplexMul(ComplexNum a, ComplexNum b);
ComplexNum ComplexDiv(ComplexNum a, ComplexNum b);
ComplexNum ComplexAdd(ComplexNum a, ComplexNum b);
ComplexNum ComplexSub(ComplexNum a, ComplexNum b);
double ComplexAbs(ComplexNum a);
double ComplexAbsSqrt(ComplexNum a);
ComplexNum ComplexConj(ComplexNum a);
ComplexNum ComplexMinus(ComplexNum a);
ComplexNum ComplexMulCon(ComplexNum a, double k);
ComplexNum ComplexNumMulExpjTheta(ComplexNum a, double Theta);
int ComplexNumMatrixMul(ComplexNum** Output, ComplexNum** a, int a_size_M, int a_size_N, int a_TranspositionEn, ComplexNum** b, int b_size_M, int b_size_N, int b_TranspositionEn);

#endif
#pragma once
