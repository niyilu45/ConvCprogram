#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fft.h"
#include "conv.h"
#include "ComplexNum.h"
int FileReadDouble(double** out, char* fileName);
int main(int argc, char* argv[])
{
    int i = 0;
    int len1, len2, len3;
    double* in1;
    double* in2;
    double* out;


    len1 = FileReadDouble(&in1, "SeqConvInA.txt");
    len2 = FileReadDouble(&in2, "SeqConvInB.txt");
    len3 = FileReadDouble(&out, "SeqConvOutC.txt");
    int convOutLen = len1 + len2 - 1;
    ComplexNum* cIn1 = (ComplexNum *)calloc(len1, sizeof(ComplexNum));
    ComplexNum* cIn2 = (ComplexNum *)calloc(len2, sizeof(ComplexNum));
    for (int i = 0; i < len1; i++) {
        cIn1[i].re = in1[i];
    }
    for (int i = 0; i < len2; i++) {
        cIn2[i].re = in2[i];
    }
    double* out1 = (double *)calloc(convOutLen, sizeof(double));
    ComplexNum* out2 = (ComplexNum *)calloc(convOutLen, sizeof(ComplexNum));
    ComplexNum* out3 = (ComplexNum *)calloc(convOutLen, sizeof(ComplexNum));

    DirectConv(out1, in1, len1, in2, len2);
    FFTConvComplex(out2, cIn1, len1, cIn2, len2);
    SegFFTConvComplex(out3, cIn1, len1, cIn2, len2);
    double maxDiff1 = 0;
    double maxDiff2 = 0;
    double maxDiff3 = 0;
    double diff;
    for (int i = 0; i < convOutLen; i++) {
        diff = abs(out[i] - out1[i]);
        if (diff > maxDiff1) {
            maxDiff1 = diff;
        }
        diff = abs(out[i] - out2[i].re);
        if (diff > maxDiff2) {
            maxDiff2 = diff;
        }
        diff = abs(out[i] - out3[i].re);
        if (diff > maxDiff3) {
            maxDiff3 = diff;
        }
    }
    printf("maxDiff1=%f\n", maxDiff1);
    printf("maxDiff2=%f\n", maxDiff2);
    printf("maxDiff3=%f\n", maxDiff3);


    free(cIn1);
    free(cIn2);
    free(out1);
    free(out2);
    free(out3);

    return 0;
}
int FileReadDouble(double** out, char* fileName){
    static FILE* Fid = NULL;
    int RawDataLen;
    double* RawDataRead;

    // 1) Calc the length of data in x.raw.
    Fid = fopen(fileName, "r");
    if(!Fid){ // Can not open the file.
        exit(0);
    }
    fseek(Fid, 0, SEEK_END); // Move to the end of file.
    RawDataLen = ftell(Fid) / sizeof(double); // ftell return size by bytes.
    RawDataLen--;
    rewind(Fid); // Move to the start of file to read data.
    if(RawDataLen <= 0){ // the raw file has no data.
        exit(0);
    }

    // 3) if FIFO has space, put data to FIFO.
    RawDataRead = (double*)malloc(sizeof(double) * RawDataLen);
    for (int i = 0; i < RawDataLen; i++) {
        fscanf(Fid, "%lf", RawDataRead+i);
    }
    *out = RawDataRead;

    return RawDataLen;
}
