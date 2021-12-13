#include "conv.h"
// This part is for Conv
void ConvDouble(double* output, double* input){
    int len1 = 11;
    int len2 = 4;
    int convOutLen = len1 + len2 - 1;
    double in1[] = {0.1419        , 0.4218 , 0.9157 , 0.7922   , 0.9595 , 0.6557
                                      , 0.0357 , 0.8491 , 0.9340   , 0.6787 , 0.7577};
    double in2[] = {0.7431        , 0.3922 , 0.6555 , 0.1712};
    ComplexNum* cin1 = (ComplexNum *)calloc(len1, sizeof(ComplexNum));
    ComplexNum* cin2 = (ComplexNum *)calloc(len2, sizeof(ComplexNum));
    ComplexNum* out = (ComplexNum *)calloc(convOutLen, sizeof(ComplexNum));
    for (int i = 0; i < len1; i++) {
        cin1[i].re = in1[i];
    }
    for (int i = 0; i < len2; i++) {
        cin2[i].re = in2[i];
    }

    /*DirectConv(out, cin1, len1, cin2, len2);*/
    /*FFTConvComplex(out, cin1, len1, cin2, len2);*/
    SegFFTConvComplex(out, cin1, len1, cin2, len2);

    for (int i = 0; i < convOutLen; i++) {
        printf("cout[%d].re = %.4f\n", i, out[i].re);
    }

    free(out);
    free(cin1);
    free(cin2);
    return;
}
void DirectConv(double* output, double* input1, int inputLen1, double* input2, int inputLen2){
    int convOutLen = inputLen1 + inputLen2 - 1;
    double* convOutTmp;
    convOutTmp = (double *)calloc(convOutLen, sizeof(double));
    for(int i=0;i<inputLen1;i++){
        for (int j = 0; j < inputLen2; j++) {
            convOutTmp[i+j] += input1[i] * input2[j];
        }
    }
    memcpy(output, convOutTmp, sizeof(double) * convOutLen);
    free(convOutTmp);
    convOutTmp = NULL;
    return;
}

void SegFFTConvComplex(ComplexNum* output, ComplexNum* input1, int inputLen1, ComplexNum* input2, int inputLen2){
    ComplexNum* in1;
    ComplexNum* in2;
    int n1, n2;
    // 1) set in1 as the longer one, in2 as the shorter one.
    if (inputLen1 > inputLen2) {
        in1 = input1;
        in2 = input2;
        n1 = inputLen1;
        n2 = inputLen2;
    } else {
        in1 = input2;
        in2 = input1;
        n1 = inputLen2;
        n2 = inputLen1;
    }

    // 2) choose the fftSize, and split the in1 by length of block.
    int fftSize;
    SegFFTConvComplexity(&fftSize, n1, n2);
    int block = fftSize + 1 - n2;
    int segs = (int)ceil((double)n1 / (double)block);

    // 3) conv.
    // 3) -0] malloc.
    int idx = 0;
    int convOutLen         = n1 + n2 - 1;
    int convOutLenTmp      = convOutLen + segs * block - n1; // this maybe different with convOutLen.
    ComplexNum* in1Tmp     = (ComplexNum *)calloc(segs * block, sizeof(ComplexNum));
    ComplexNum* X1         = (ComplexNum *)malloc(fftSize * sizeof(ComplexNum));
    ComplexNum* X2         = (ComplexNum *)calloc(fftSize, sizeof(ComplexNum));
    ComplexNum* Y          = (ComplexNum *)malloc(fftSize * sizeof(ComplexNum));
    ComplexNum* y          = (ComplexNum *)malloc(fftSize * sizeof(ComplexNum));
    ComplexNum* convOutTmp = (ComplexNum *)calloc(convOutLenTmp, sizeof(ComplexNum));

    // this can be considered usede parallel for.
    memcpy(in1Tmp, in1, n1 * sizeof(ComplexNum));
    memcpy(X2, in2, n2 * sizeof(ComplexNum));
    FFT(X2, X2, fftSize);
    for (int segCnt = 0; segCnt < segs; segCnt++) {
        // 3) -1] Input fft.
        memset(X1, 0, fftSize * sizeof(ComplexNum));
        memcpy(X1, in1Tmp + idx, block * sizeof(ComplexNum));
        FFT(X1, X1, fftSize);

        // 3) -2] Frequency domain multiplication.
        double mulTmpA, mulTmpB, mulTmpC;
        for (int i = 0; i < fftSize; i++) {
            mulTmpA = (X1[i].re + X1[i].im) * X2[i].re;
            mulTmpB = (X2[i].re + X2[i].im) * X1[i].im;
            mulTmpC = (X1[i].im - X1[i].re) * X2[i].im;
            Y[i].re = mulTmpA - mulTmpB;
            Y[i].im = mulTmpB - mulTmpC;
        }

        // 3) -3] output ifft.
        IFFT(y, Y, fftSize);

        // 3) -4] Combine.
        ComplexNum* ptr1 = convOutTmp + idx;
        ComplexNum* ptr2 = y;
        for (int i = 0; i < fftSize; i++) {
            (*ptr1).re += (*ptr2).re;
            (*ptr1).im += (*ptr2).im;
            ptr1++;
            ptr2++;
        }
        idx += block;
    }
    memcpy(output, convOutTmp, sizeof(ComplexNum) * convOutLen);

    free(in1Tmp);
    free(X1);
    free(X2);
    free(Y);
    free(y);
    free(convOutTmp);

    return;
}

enum convMethod FindBestConvMethod(int n1, int n2){
    int c1 = DirectConvComplexity(n1, n2);
    int bestFFTSize;
    int c2 = SegFFTConvComplexity(&bestFFTSize, n1, n2);
    return (c2 < c1)?segFFTConv:directConv;
}
int DirectConvComplexity(int n1, int n2){
    return n1 * n2;
}
int SegFFTConvComplexity(int* bestFFTSize, int n1, int n2){
    if (n1 < n2) {
        int tmp = n1;
        n1 = n2;
        n2 = tmp;
    }
    int minFFTOrder = CalcFFTOrder(1+n2-1);
    int maxFFTOrder = CalcFFTOrder(n1+n2-1);
    int c = INT_MAX; // times of multiply
    *bestFFTSize = CalcFFTLen(1+n2-1);
    int fftOrder = minFFTOrder;
    while (fftOrder <= maxFFTOrder) {
        // 1) calc depart of n1.
        int fftSize = 1 << fftOrder;
        int block = fftSize + 1 - n2; // n1 was depart by length of block.
        int segs = ceil(n1 / block); // n1 was depart into n(segs) blocks.
        // 2) find the min complexity.
        int doubleMulPerComplexMul = 3;
        int totalFFTTimes = 2 * segs + 1;
        int cPerFFT = (fftSize * fftOrder) >> 1;
        int complexMulTimes = segs * fftSize;
        int complexity = doubleMulPerComplexMul * (totalFFTTimes * cPerFFT + complexMulTimes);
        if (complexity < c) {
            c = complexity;
            *bestFFTSize = fftSize;
        }
        fftOrder++;
    }
    return c;
}
