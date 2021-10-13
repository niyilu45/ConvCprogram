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
    int i,j;
    int convOutLen = inputLen1 + inputLen2 - 1;
    double* convOutTmp;
    convOutTmp = (double *)calloc(convOutLen, sizeof(double));
    for(i=0;i<inputLen1;i++){
        for (j = 0; j < inputLen2; j++) {
            convOutTmp[i+j] += input1[i] * input2[j];
        }
    }
    memcpy(output, convOutTmp, sizeof(double) * convOutLen);
    free(convOutTmp);
    convOutTmp = NULL;
    return;
}
void FFTConvComplex(ComplexNum* output, ComplexNum* input1, int inputLen1, ComplexNum* input2, int inputLen2){
    int convOutLen = inputLen1 + inputLen2 - 1;
    int fftLen = CalcFFTSize(convOutLen);
    ComplexNum* x1 = (ComplexNum *)calloc(fftLen, sizeof(ComplexNum));
    ComplexNum* x2 = (ComplexNum *)calloc(fftLen, sizeof(ComplexNum));
    memcpy(x1, input1, sizeof(ComplexNum) * inputLen1);
    memcpy(x2, input2, sizeof(ComplexNum) * inputLen2);
    ComplexNum* X1 = FFT(x1, fftLen);
    ComplexNum* X2 = FFT(x2, fftLen);
    ComplexNum* Y = (ComplexNum *)calloc(fftLen, sizeof(ComplexNum));
    for (int i = 0; i < fftLen; i++) {
        double A = (X1[i].re + X1[i].im) * X2[i].re;
        double B = (X2[i].re + X2[i].im) * X1[i].im;
        double C = (X1[i].im - X1[i].re) * X2[i].im;
        Y[i].re = A - B;
        Y[i].im = B - C;
    }
    ComplexNum* y = IFFT(Y, fftLen);
    memcpy(output, y, sizeof(ComplexNum) * convOutLen);

    free(x1);
    free(x2);
    free(X1);
    free(X2);
    free(Y);
    free(y);
    return;
}
void SegFFTConvComplex(ComplexNum* output, ComplexNum* input1, int inputLen1, ComplexNum* input2, int inputLen2){
    ComplexNum* in1;
    ComplexNum* in2;
    int n1, n2;
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
    int convOutLen = n1 + n2 - 1;
    int fftLen = CalcFFTSize(2 * n2 - 1);
    int block = fftLen + 1 - n2;
    int segs = (int)ceil((double)n1 / (double)block);

    ComplexNum* x1 = (ComplexNum *)calloc(fftLen, sizeof(ComplexNum));
    ComplexNum* x2 = (ComplexNum *)calloc(fftLen, sizeof(ComplexNum));
    memcpy(x2, input2, sizeof(ComplexNum) * inputLen2);
    ComplexNum* X2 = FFT(x2, fftLen);
    ComplexNum* Y = (ComplexNum *)calloc(fftLen, sizeof(ComplexNum));
    int outTmpLen = segs*block+n2-1;
    ComplexNum* outTmp = (ComplexNum *)calloc(outTmpLen, sizeof(ComplexNum));
    ComplexNum* outTmpPtr = outTmp;
    int curIn1Idx = 0;
    for (int i = 0; i < segs; i++) {
        if (curIn1Idx + block < n1) {
            memcpy(x1, in1 + curIn1Idx, sizeof(ComplexNum) * block);
            memset(x1 + block, 0, sizeof(ComplexNum) * (fftLen - block));
        }
        else{
            memcpy(x1, in1 + curIn1Idx, sizeof(ComplexNum) * (n1 - curIn1Idx));
            memset(x1 + n1 - curIn1Idx, 0, sizeof(ComplexNum) * (fftLen - n1 + curIn1Idx));
        }
        ComplexNum* X1 = FFT(x1, fftLen);
        for (int i = 0; i < fftLen; i++) {
            double A = (X1[i].re + X1[i].im) * X2[i].re;
            double B = (X2[i].re + X2[i].im) * X1[i].im;
            double C = (X1[i].im - X1[i].re) * X2[i].im;
            Y[i].re = A - B;
            Y[i].im = B - C;
        }
        ComplexNum* y = IFFT(Y ,fftLen);
        ComplexNum* yPtr = y;
        for (int i = 0; i < fftLen-block; i++) {
            outTmpPtr[i].re += y[i].re;
            outTmpPtr[i].im += y[i].im;
        }
        memcpy(&outTmpPtr[fftLen-block], &y[fftLen-block], sizeof(ComplexNum) * block);
        outTmpPtr += block;
        curIn1Idx += block;
        free(X1);
        free(y);
    }

    memcpy(output, outTmp, sizeof(ComplexNum) * convOutLen);
    free(x1);
    free(x2);
    free(X2);
    free(Y);
    free(outTmp);
}
int CalcFFTSize(int InputLen){
    int FFTLen = 1 << (int)ceil((log(InputLen)/log(2)));
    return FFTLen;
}
