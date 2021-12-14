#ifndef __CONV_H
#define __CONV_H
#include "ComplexNum.h"
#include "fft.h"
enum convMethod {
    directConv = 0,
    segFFTConv
};
void ConvDouble(double* output, double* input1, int inputLen1, double* input2, int inputLen2);
enum convMethod FindBestConvMethod(int n1, int n2);
void DirectConv(double* output, double* input1, int inputLen1, double* input2, int inputLen2);
void SegFFTConvComplex(ComplexNum* output, ComplexNum* input1, int inputLen1, ComplexNum* input2, int inputLen2);
int DirectConvComplexity(int n1, int n2);
int SegFFTConvComplexity(int* bestFFTSize, int n1, int n2);
#endif
