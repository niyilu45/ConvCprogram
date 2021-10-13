#include "ComplexNum.h"

ComplexNum ComplexNumReset(void)
{
    ComplexNum Z;
    Z.re = 0;
    Z.im = 0;

    return Z;
}
ComplexNum ComplexMul(ComplexNum a, ComplexNum b)
{
    ComplexNum Mul;
    Mul.re = a.re * b.re - a.im * b.im;
    Mul.im = a.re * b.im + a.im * b.re;
    return Mul;
}
ComplexNum ComplexDiv(ComplexNum a, ComplexNum b)
{
    double delta = b.re * b.re + b.im * b.im;
    ComplexNum Div;
    Div.re = (a.re * b.re + a.im * b.im) / delta;
    Div.im = (b.re * a.im - a.re * b.im) / delta;
    return Div;
}
ComplexNum ComplexAdd(ComplexNum a, ComplexNum b)
{
    ComplexNum Add;
    Add.re = a.re + b.re;
    Add.im = a.im + b.im;
    return Add;
}
ComplexNum ComplexSub(ComplexNum a, ComplexNum b)
{
    ComplexNum Sub;
    Sub.re = a.re - b.re;
    Sub.im = a.im - b.im;
    return Sub;
}
double ComplexAbs(ComplexNum a)
{
    double abs = 0;
    abs = a.re * a.re + a.im * a.im;
    return abs;
}
double ComplexAbsSqrt(ComplexNum a)
{
    double abssqrt = 0;
    abssqrt = sqrt(a.re * a.re + a.im * a.im);
    return abssqrt;
}
ComplexNum ComplexConj(ComplexNum a)
{
    ComplexNum Conj;
    Conj.re = a.re;
    Conj.im = -a.im;
    return Conj;
}
ComplexNum ComplexMinus(ComplexNum a)
{
    ComplexNum Minus;
    Minus.re = -a.re;
    Minus.im = -a.im;
    return Minus;
}
ComplexNum ComplexMulCon(ComplexNum a, double k)
{
    ComplexNum Mulc;
    Mulc.re = a.re * k;
    Mulc.im = a.im * k;

    return Mulc;
}
ComplexNum ComplexNumMulExpjTheta(ComplexNum a, double Theta)
{
    ComplexNum MuljT;
    double CosTmp = cos(Theta);
    double SinTmp = sin(Theta);
    MuljT.re = a.re * CosTmp - a.im * SinTmp;
    MuljT.im = a.im * CosTmp + a.re * SinTmp;

    return MuljT;
}
int ComplexNumMatrixMul(ComplexNum** Output, ComplexNum** a, int a_size_M, int a_size_N, int a_TranspositionEn, ComplexNum** b, int b_size_M, int b_size_N, int b_TranspositionEn){
    int i,j;
    int Temp;
    int a_size_M_T, a_size_N_T;
    int b_size_M_T, b_size_N_T;
    ComplexNum** a_T;
    ComplexNum** b_T;

    if(a_TranspositionEn == 1){
        a_size_M_T = a_size_N;
        a_size_N_T = a_size_M;
    }
    else{
        a_size_M_T = a_size_M;
        a_size_N_T = a_size_N;
    }
    if(b_TranspositionEn == 1){
        b_size_M_T = b_size_N;
        b_size_N_T = b_size_M;
    }
    else{
        b_size_M_T = b_size_M;
        b_size_N_T = b_size_N;
    }

    if(a_size_N_T != b_size_M_T){
        return -1;
    }

    a_T = (ComplexNum**)malloc(sizeof(ComplexNum*) * a_size_M_T);
    for(i=0;i<a_size_M_T;i++){
        a_T[i] = (ComplexNum*)malloc(sizeof(ComplexNum) * a_size_N_T);
    }
    b_T = (ComplexNum**)malloc(sizeof(ComplexNum*) * b_size_M_T);
    for(i=0;i<b_size_M_T;i++){
        b_T[i] = (ComplexNum*)malloc(sizeof(ComplexNum) * b_size_N_T);
    }

    for(i=0;i<a_size_M_T;i++){
        for(j=0;j<a_size_N_T;j++){
            a_T[i][j] = (a_TranspositionEn)?a[j][i]:a[i][j];
        }
    }
    for(i=0;i<b_size_M_T;i++){
        for(j=0;j<b_size_N_T;j++){
            b_T[i][j] = (b_TranspositionEn)?b[j][i]:b[i][j];
        }
    }

    int m,n;
    double Re,Im;
    for(m=0;m<a_size_M_T;m++){
        for(n=0;n<b_size_N_T;n++){
            Re = 0;
            Im = 0;
            for(i=0;i<a_size_N_T;i++){
                Re += a_T[m][i].re * b_T[i][n].re - a_T[m][i].im * b_T[i][n].im;
                Im += a_T[m][i].re * b_T[i][n].im + a_T[m][i].im * b_T[i][n].re;
            }
            Output[m][n].re = Re;
            Output[m][n].im = Im;
        }
    }


    for(i=0;i<a_size_M_T;i++){
        free(a_T[i]);
        a_T[i] = NULL;
    }
    free(a_T);
    a_T = NULL;
    for(i=0;i<b_size_M_T;i++){
        free(b_T[i]);
        b_T[i] = NULL;
    }
    free(b_T);
    b_T = NULL;

    return 0;
}
