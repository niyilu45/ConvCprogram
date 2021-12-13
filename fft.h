#ifndef FFT_H
#define FFT_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "ComplexNum.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif

/*
typedef struct ComplexNum{
    double re;
    double im;
}ComplexNum;
*/

typedef struct OmegaLibNode* OmegaLibList;
typedef struct OmegaLibNode{
    OmegaLibList next;
    int len;
    ComplexNum* omega;
}OmegaLibNode;


int CalcFFTLen(int InputLen);
int CalcFFTOrder(int InputLen);
void ClearOmegaLib(void); // After FFT, should clear lib.
void FFT(ComplexNum* output, ComplexNum* input, int inputLen);
void IFFT(ComplexNum* output, ComplexNum* input, int inputLen);


static void FFTDataFlip(ComplexNum* input, int inputLen);
static int CalcFFTStateNum(int inputLen);
static unsigned int BitFlipInt(register unsigned int val, int BitWidth);

// Omega Lib
static void InsertOmegaLib(OmegaLibList L, int fftLen);
static OmegaLibList FindOmegaLib(int fftLen, OmegaLibList head);
static ComplexNum CalcOmega(int i, int N);
static void CalcOmegaLib(ComplexNum* omega, int len);
static OmegaLibList DestroyOmegaLib(OmegaLibList L);

#endif
