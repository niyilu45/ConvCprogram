#ifndef __CONV_H
#define __CONV_H
#include "ComplexNum.h"
#include "fft.h"
void ConvDouble(double* output, double* input);
void DirectConv(double* output, double* input1, int inputLen1, double* input2, int inputLen2);
void FFTConvComplex(ComplexNum* output, ComplexNum* input1, int inputLen1, ComplexNum* input2, int inputLen2);
void SegFFTConvComplex(ComplexNum* output, ComplexNum* input1, int inputLen1, ComplexNum* input2, int inputLen2);
int CalcFFTSize(int InputLen);
#endif
