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
                Temp = ComplexMulInFFT(Output[j+k+HalfBlkSize], Omega[k*(FFTLen/BlkSize)]);
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

    double tA = (A.re + A.im) * B.re;
    double tB = (B.re + B.im) * A.im;
    double tC = (A.im - A.re) * B.im;
    Output.re = tA - tB;
    Output.im = tB - tC;

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

int CalFFTLen(int InputLen){
    int FFTLen = 1 << (int)ceil((log(InputLen)/log(2)));
    return FFTLen;
}
