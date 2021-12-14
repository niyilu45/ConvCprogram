#include "fft.h"

OmegaLibList omegaLib = NULL;

void FFT(ComplexNum* output, ComplexNum* input, int inputLen){
    if(inputLen <= 0){
        return;
    }

    // Check the fftLen.
    int stateNum = CalcFFTStateNum(inputLen);
    int fftLen = CalcFFTLen(inputLen);

    if(inputLen != fftLen){
        printf("\r\n[Error] in %s! Do not support none pow of 2 points!", __func__);
        exit(0);
    }

    // 1) Prepare omega.
    ComplexNum* omega;
    OmegaLibList curOmegaLib;

    omegaLib = InsertOmegaLib(omegaLib, fftLen);
    curOmegaLib = FindOmegaLib(fftLen, omegaLib);
    omega = curOmegaLib->omega;

    // 2) Change the input Order.
    ComplexNum* outputTmp = (ComplexNum *)malloc(fftLen * sizeof(ComplexNum));
    memcpy(outputTmp, input, fftLen * sizeof(ComplexNum));
    FFTDataFlip(outputTmp, fftLen);

    // 3) Main Part.
    int blkSize, halfBlkSize;
    ComplexNum inA, inB;
    ComplexNum mulTmp;
    ComplexNum addTmp, subTmp;
    double mulTmpA, mulTmpB, mulTmpC;
    for(int i=0;i<stateNum;i++){
        blkSize = 2<<i;
        halfBlkSize = 1<<i;
        for(int j=0;j<fftLen;j+=blkSize){
            for(int k=0;k<halfBlkSize;k++){
                // complex multiply
                inA = outputTmp[j+k+halfBlkSize];
                inB = omega[k*(fftLen/blkSize)];
                mulTmpA = (inA.re + inA.im) * inB.re;
                mulTmpB = (inB.re + inB.im) * inA.im;
                mulTmpC = (inA.im - inA.re) * inB.im;
                mulTmp.re = mulTmpA - mulTmpB;
                mulTmp.im = mulTmpB - mulTmpC;

                // complex add and sub.
                addTmp.re = outputTmp[j+k].re+mulTmp.re;
                addTmp.im = outputTmp[j+k].im+mulTmp.im;
                subTmp.re = outputTmp[j+k].re-mulTmp.re;
                subTmp.im = outputTmp[j+k].im-mulTmp.im;
                outputTmp[j+k] = addTmp;
                outputTmp[j+k+halfBlkSize] = subTmp;
            }
        }

    }

    memcpy(output, outputTmp, fftLen * sizeof(ComplexNum));
    free(outputTmp);;
    return;
}

void IFFT(ComplexNum* output, ComplexNum* input, int inputLen){
    if(inputLen <= 0){
        return;
    }
    // Check the fftLen.
    int stateNum = CalcFFTStateNum(inputLen);
    int fftLen = CalcFFTLen(inputLen);

    if(inputLen != fftLen){
        printf("\r\nFFT Err! Do not support none pow of 2 points!");
        exit(0);
    }

    ComplexNum* inputTmp;
    inputTmp = (ComplexNum*)malloc(sizeof(ComplexNum) * fftLen);
    for(int i=0;i<fftLen;i++){
        inputTmp[i].re = input[i].re;
        inputTmp[i].im = -1*input[i].im;
    }

    FFT(inputTmp, inputTmp, fftLen);
    memcpy(output, inputTmp, fftLen * sizeof(ComplexNum));

    for(int i=0;i<fftLen;i++){
        output[i].re /= fftLen;
        output[i].im /= (-fftLen);
    }

    free(inputTmp);
    inputTmp = NULL;

    return;
}


static ComplexNum CalcOmega(int i, int N){
    ComplexNum omega;
    omega.re = cos(2*PI*i/N);
    omega.im = -sin(2*PI*i/N);

    return omega;
}

static void FFTDataFlip(ComplexNum* input, int inputLen){
    int BitWidth = CalcFFTStateNum(inputLen);
    unsigned int idx;

    ComplexNum* outputTmp = (ComplexNum *)malloc(inputLen * sizeof(ComplexNum));
    
    for(int i=0;i<inputLen;i++){
        idx = BitFlipInt(i, BitWidth);
        outputTmp[idx] = input[i];
    }

    memcpy(input, outputTmp, inputLen * sizeof(ComplexNum));
    free(outputTmp);
    return;
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
        printf("[Error] in %s! Do not support this yet!\n", __func__);
        exit(0);
    }

    return val;
}

static OmegaLibList FindOmegaLib(int fftLen, OmegaLibList head){
    OmegaLibList L = head;
    while(L != NULL){
        if(L->len == fftLen){
            return L;
        }
        L = L->next;
    }

    return NULL;
}

static OmegaLibList InsertOmegaLib(OmegaLibList L, int fftLen){
    // 1) if find the end of the list, create new.
    if(L == NULL){
        L = (OmegaLibList)calloc(1, sizeof(OmegaLibNode));
        L->next = NULL;
        L->len = fftLen;
        L->omega = (ComplexNum *)malloc(fftLen * sizeof(ComplexNum));
        CalcOmegaLib(L->omega, fftLen);
        return L;
    }

    // 2) if this len is already exist, skip.
    if(L->len == fftLen){
        return L;
    }

    // 3) find next list node.
    L->next = InsertOmegaLib(L->next, fftLen);
    return L;
}

static void CalcOmegaLib(ComplexNum* omega, int len){
    for(int i=0;i<len;i++){
        omega[i].re = cos(2*PI*i/len);
        omega[i].im = sin(-2*PI*i/len);
    }
    return;
}

static OmegaLibList DestroyOmegaLib(OmegaLibList L){
    if(L == NULL){
        return NULL;
    }

    L->next = DestroyOmegaLib(L->next);
    free(L->omega);
    L->omega = NULL;
    free(L);

    return NULL;
}

void ClearOmegaLib(void){
    omegaLib = DestroyOmegaLib(omegaLib);
    return;
}

static int CalcFFTStateNum(int inputLen){
    int stateNum = 0;
    inputLen--;
    while (inputLen) {
        inputLen >>= 1;
        stateNum++;
    }
    return stateNum;
}

int CalcFFTOrder(int inputLen){
    int fftOrder = 0;
    inputLen--;
    while (inputLen) {
        inputLen >>= 1;
        fftOrder++;
    }
    return fftOrder;
}

int CalcFFTLen(int inputLen){
    int fftOrder = 0;
    inputLen--;
    while (inputLen) {
        inputLen >>= 1;
        fftOrder++;
    }
    return 1 << fftOrder;
}
