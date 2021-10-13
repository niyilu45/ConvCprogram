#include "fft.h"

OmegaLibList OmegaLib = NULL;

ComplexNum* FFT(ComplexNum* Input, int InputLen){
    if(InputLen <= 0){
        return NULL;
    }

    int i,j,k;
    // Check the FFTLen.
    int StateNum = CalFFTStateNum(InputLen);
    int FFTLen = CalFFTLen(InputLen);

    if(InputLen != FFTLen){
        printf("\r\nFFT Err! Do not support none pow of 2 points!");
        exit(0);
    }

    // 1) Prepare Omega.
    ComplexNum* Omega;
    OmegaLibList CurOmegaLib;

    while((CurOmegaLib = FindOmegaLib(FFTLen, OmegaLib)) == NULL){
        OmegaLib = InsertOmegaLib(FFTLen, OmegaLib);
    }
    Omega = CurOmegaLib->Omega;

    // 2) Change the Input Order.
    ComplexNum* Output;
    Output = FFTDataFlip(Input, FFTLen);

    // 3) Main Part.
    int BlkSize, HalfBlkSize;
    ComplexNum Temp;
    ComplexNum AddTemp, SubTemp;
    for(i=0;i<StateNum;i++){
        BlkSize = 2<<i;
        HalfBlkSize = 1<<i;
        for(j=0;j<FFTLen;j+=BlkSize){
            for(k=0;k<HalfBlkSize;k++){
                Temp = ComplexMulInFFT(Output[j+k+HalfBlkSize], Omega[k*FFTLen/BlkSize]);
                //AddTemp = ComplexAddInFFT(Output[j+k],Temp);
                //SubTemp = ComplexSubInFFT(Output[j+k],Temp);
                AddTemp.re = Output[j+k].re+Temp.re;
                AddTemp.im = Output[j+k].im+Temp.im;
                SubTemp.re = Output[j+k].re-Temp.re;
                SubTemp.im = Output[j+k].im-Temp.im;
                Output[j+k] = AddTemp;
                Output[j+k+HalfBlkSize] = SubTemp;
            }
        }

    }

    return Output;
}

ComplexNum* IFFT(ComplexNum* Input, int InputLen){
    if(InputLen <= 0){
        return NULL;
    }
    int i;
    // Check the FFTLen.
    int StateNum = CalFFTStateNum(InputLen);
    int FFTLen = CalFFTLen(InputLen);

    if(InputLen != FFTLen){
        printf("\r\nFFT Err! Do not support none pow of 2 points!");
        exit(0);
    }

    ComplexNum* InputTemp;
    InputTemp = (ComplexNum*)malloc(sizeof(ComplexNum) * FFTLen);
    for(i=0;i<FFTLen;i++){
        InputTemp[i].re = Input[i].re;
        InputTemp[i].im = -1*Input[i].im;
    }

    ComplexNum* Output;
    Output = FFT(InputTemp, FFTLen);

    for(i=0;i<FFTLen;i++){
        Output[i].re /= FFTLen;
        Output[i].im /= (-FFTLen);
    }

    free(InputTemp);
    InputTemp = NULL;

    return Output;
}


static ComplexNum CalOmega(int i, int N){
    ComplexNum Omega;
    Omega.re = cos(2*PI*i/N);
    Omega.im = -sin(2*PI*i/N);

    return Omega;
}

static ComplexNum ComplexAddInFFT(ComplexNum A, ComplexNum B){
    ComplexNum Output;

    Output.re = A.re + B.re;
    Output.im = A.im + B.im;

    return Output;
}
static ComplexNum ComplexSubInFFT(ComplexNum A, ComplexNum B){
    ComplexNum Output;

    Output.re = A.re - B.re;
    Output.im = A.im - B.im;

    return Output;
}

static ComplexNum ComplexMulInFFT(ComplexNum A, ComplexNum B){
    ComplexNum Output;

    Output.re = A.re * B.re - A.im * B.im;
    Output.im = A.re * B.im + A.im * B.re;

    return Output;
}

static ComplexNum* FFTDataFlip(ComplexNum* Input, int InputLen){
    int i;
    ComplexNum* FlipedData;
    int BitWidth = CalFFTStateNum(InputLen);
    unsigned int Idx;

    FlipedData = (ComplexNum*)malloc(sizeof(ComplexNum) * InputLen);
    for(i=0;i<InputLen;i++){
        Idx = BitFlipInt(i, BitWidth);
        FlipedData[Idx] = Input[i];
    }

    return FlipedData;
}

static int CalFFTStateNum(int InputLen){
    int StateNum = (int)ceil((log(InputLen)/log(2)));
    return StateNum;
}

static unsigned int BitFlipInt(register unsigned int val, int BitWidth){
    if(sizeof(unsigned int) == 4){
        val = (val & 0xaaaaaaaa) >> 1 | (val & 0x55555555) << 1;
        val = (val & 0xcccccccc) >> 2 | (val & 0x33333333) << 2;
        val = (val & 0xf0f0f0f0) >> 4 | (val & 0x0f0f0f0f) << 4;
        val = (val & 0xff00ff00) >> 8 | (val & 0x00ff00ff) << 8;
        val = val >> 16 | val << 16;
        val >>= (32-BitWidth);
    }
    else{
        printf("BitFlip Err! Do not support this yet!\n");
        exit(0);
    }

    return val;
}

static OmegaLibList FindOmegaLib(int FFTLen, OmegaLibList Head){
    OmegaLibList L;
    L = Head;
    while(L != NULL){
        if(L->Len == FFTLen){
            return L;
        }
        L = L->Next;
    }

    return NULL;
}

static OmegaLibList InsertOmegaLib(int FFTLen, OmegaLibList L){
    if(L == NULL){
        L = (OmegaLibList)calloc(1, sizeof(OmegaLibNode));
        L->Next = NULL;
        L->Len = FFTLen;
        L->Omega = CalcOmegaLib(FFTLen);
        return L;
    }

    if(L->Len == FFTLen){
        return L;
    }

    L->Next = InsertOmegaLib(FFTLen, L->Next);

    return L;
}

static ComplexNum* CalcOmegaLib(int FFTLen){
    int i;
    ComplexNum* Omega;
    Omega = (ComplexNum*)malloc(sizeof(ComplexNum) * FFTLen);

    for(i=0;i<FFTLen;i++){
        Omega[i].re = cos(2*PI*i/FFTLen);
        Omega[i].im = sin(-2*PI*i/FFTLen);
    }

    return Omega;
}

void ClearOmegaLib(void){
    OmegaLib = DestroyOmegaLib(OmegaLib);
    return;
}

static OmegaLibList DestroyOmegaLib(OmegaLibList L){
    if(L == NULL){
        return NULL;
    }

    L->Next = DestroyOmegaLib(L->Next);
    free(L->Omega);
    L->Omega = NULL;
    free(L);

    return NULL;
}

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

    [>DirectConv(out, cin1, len1, cin2, len2);<]
    [>FFTConvComplex(out, cin1, len1, cin2, len2);<]
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
    int fftLen = CalFFTLen(convOutLen);
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
    int fftLen = CalFFTLen(2 * n2 - 1);
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
            [>memset(x1 + block, 0, sizeof(ComplexNum) * (fftLen - block)); // This is not needed<]
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
int CalFFTLen(int InputLen){
    int FFTLen = 1 << (int)ceil((log(InputLen)/log(2)));
    return FFTLen;
}
