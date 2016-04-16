/**********************************************************************
 
 fft.cpp
 
 
 This class is a C++ wrapper for original code
 written by Dominic Mazzoni in September 2000
 
 This file contains a few FFT routines, including a real-FFT
 routine that is almost twice as fast as a normal complex FFT,
 and a power spectrum routine which is more convenient when
 you know you don't care about phase information.  It now also
 contains a few basic windowing functions.
 
 Some of this code was based on a free implementation of an FFT
 by Don Cross, available on the web at:
 
 http://www.intersrv.com/~dcross/fft.html
 
 The basic algorithm for his code was based on Numerical Recipes
 in Fortran.  I optimized his code further by reducing array
 accesses, caching the bit reversal table, and eliminating
 float-to-double conversions, and I added the routines to
 calculate a real FFT and a real power spectrum.
 
 Note: all of these routines use single-precision floats.
 I have found that in practice, floats work well until you
 get above 8192 samples.  If you need to do a larger FFT,
 you need to use doubles.
 
 **********************************************************************/

#include "fft.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int **gFFTBitTable = NULL;
const int MaxFastBits = 16;

int IsPowerOfTwo(int x)
{
    if (x < 2)
        return false;
    
    if (x & (x - 1))             /* Thanks to 'byang' for this cute trick! */
        return false;
    
    return true;
}

int NumberOfBitsNeeded(int PowerOfTwo)
{
    int i;
    
    if (PowerOfTwo < 2) {
        fprintf(stderr, "Error: FFT called with size %d\n", PowerOfTwo);
        exit(1);
    }
    
    for (i = 0;; i++)
        if (PowerOfTwo & (1 << i))
            return i;
}

int ReverseBits(int index, int NumBits)
{
    int i, rev;
    
    for (i = rev = 0; i < NumBits; i++) {
        rev = (rev << 1) | (index & 1);
        index >>= 1;
    }
    
    return rev;
}

void InitFFT()
{
    gFFTBitTable = new int *[MaxFastBits];
    
    int len = 2;
    for (int b = 1; b <= MaxFastBits; b++) {
        
        gFFTBitTable[b - 1] = new int[len];
        
        for (int i = 0; i < len; i++)
            gFFTBitTable[b - 1][i] = ReverseBits(i, b);
        
        len <<= 1;
    }
}

inline int FastReverseBits(int i, int NumBits)
{
    if (NumBits <= MaxFastBits)
        return gFFTBitTable[NumBits - 1][i];
    else
        return ReverseBits(i, NumBits);
}

/*
 * Complex Fast Fourier Transform
 */

void FFT(int NumSamples,
         bool InverseTransform,
         float *RealIn, float *ImagIn, float *RealOut, float *ImagOut)
{
    int NumBits;                 /* Number of bits needed to store indices */
    int i, j, k, n;
    int BlockSize, BlockEnd;
    
    double angle_numerator = 2.0 * M_PI;
    float tr, ti;                /* temp real, temp imaginary */
    
    if (!IsPowerOfTwo(NumSamples)) {
        fprintf(stderr, "%d is not a power of two\n", NumSamples);
        exit(1);
    }
    
    if (!gFFTBitTable)
        InitFFT();
    
    if (InverseTransform)
        angle_numerator = -angle_numerator;
    
    NumBits = NumberOfBitsNeeded(NumSamples);
    
    /*
     **   Do simultaneous data copy and bit-reversal ordering into outputs...
     */
    
    for (i = 0; i < NumSamples; i++) {
        j = FastReverseBits(i, NumBits);
        RealOut[j] = RealIn[i];
        ImagOut[j] = (ImagIn == NULL) ? 0.0 : ImagIn[i];
    }
    
    /*
     **   Do the FFT itself...
     */
    
    BlockEnd = 1;
    for (BlockSize = 2; BlockSize <= NumSamples; BlockSize <<= 1) {
        
        double delta_angle = angle_numerator / (double) BlockSize;
        
        float sm2 = sin(-2 * delta_angle);
        float sm1 = sin(-delta_angle);
        float cm2 = cos(-2 * delta_angle);
        float cm1 = cos(-delta_angle);
        float w = 2 * cm1;
        float ar0, ar1, ar2, ai0, ai1, ai2;
        
        for (i = 0; i < NumSamples; i += BlockSize) {
            ar2 = cm2;
            ar1 = cm1;
            
            ai2 = sm2;
            ai1 = sm1;
            
            for (j = i, n = 0; n < BlockEnd; j++, n++) {
                ar0 = w * ar1 - ar2;
                ar2 = ar1;
                ar1 = ar0;
                
                ai0 = w * ai1 - ai2;
                ai2 = ai1;
                ai1 = ai0;
                
                k = j + BlockEnd;
                tr = ar0 * RealOut[k] - ai0 * ImagOut[k];
                ti = ar0 * ImagOut[k] + ai0 * RealOut[k];
                
                RealOut[k] = RealOut[j] - tr;
                ImagOut[k] = ImagOut[j] - ti;
                
                RealOut[j] += tr;
                ImagOut[j] += ti;
            }
        }
        
        BlockEnd = BlockSize;
    }
    
    /*
     **   Need to normalize if inverse transform...
     */
    
    if (InverseTransform) {
        float denom = (float) NumSamples;
        
        for (i = 0; i < NumSamples; i++) {
            RealOut[i] /= denom;
            ImagOut[i] /= denom;
        }
    }
}

/*
 * Real Fast Fourier Transform
 *
 * This function was based on the code in Numerical Recipes in C.
 * In Num. Rec., the inner loop is based on a single 1-based array
 * of interleaved real and imaginary numbers.  Because we have two
 * separate zero-based arrays, our indices are quite different.
 * Here is the correspondence between Num. Rec. indices and our indices:
 *
 * i1  <->  real[i]
 * i2  <->  imag[i]
 * i3  <->  real[n/2-i]
 * i4  <->  imag[n/2-i]
 */

void RealFFT(int NumSamples, float *RealIn, float *RealOut, float *ImagOut)
{
    int Half = NumSamples / 2;
    int i;
    
    float theta = M_PI / Half;
    
    float *tmpReal = new float[Half];
    float *tmpImag = new float[Half];
    
    for (i = 0; i < Half; i++) {
        tmpReal[i] = RealIn[2 * i];
        tmpImag[i] = RealIn[2 * i + 1];
    }
    
    FFT(Half, 0, tmpReal, tmpImag, RealOut, ImagOut);
    
    float wtemp = float (sin(0.5 * theta));
    
    float wpr = -2.0 * wtemp * wtemp;
    float wpi = float (sin(theta));
    float wr = 1.0 + wpr;
    float wi = wpi;
    
    int i3;
    
    float h1r, h1i, h2r, h2i;
    
    for (i = 1; i < Half / 2; i++) {
        
        i3 = Half - i;
        
        h1r = 0.5 * (RealOut[i] + RealOut[i3]);
        h1i = 0.5 * (ImagOut[i] - ImagOut[i3]);
        h2r = 0.5 * (ImagOut[i] + ImagOut[i3]);
        h2i = -0.5 * (RealOut[i] - RealOut[i3]);
        
        RealOut[i] = h1r + wr * h2r - wi * h2i;
        ImagOut[i] = h1i + wr * h2i + wi * h2r;
        RealOut[i3] = h1r - wr * h2r + wi * h2i;
        ImagOut[i3] = -h1i + wr * h2i + wi * h2r;
        
        wr = (wtemp = wr) * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
    }
    
    RealOut[0] = (h1r = RealOut[0]) + ImagOut[0];
    ImagOut[0] = h1r - ImagOut[0];
    
    delete[]tmpReal;
    delete[]tmpImag;
}

/*
 * PowerSpectrum
 *
 * This function computes the same as RealFFT, above, but
 * adds the squares of the real and imaginary part of each
 * coefficient, extracting the power and throwing away the
 * phase.
 *
 * For speed, it does not call RealFFT, but duplicates some
 * of its code.
 */

void PowerSpectrum(int NumSamples, float *In, float *Out)
{
    int Half = NumSamples / 2;
    int i;
    
    float theta = M_PI / Half;
    
    float *tmpReal = new float[Half];
    float *tmpImag = new float[Half];
    float *RealOut = new float[Half];
    float *ImagOut = new float[Half];
    
    for (i = 0; i < Half; i++) {
        tmpReal[i] = In[2 * i];
        tmpImag[i] = In[2 * i + 1];
    }
    
    FFT(Half, 0, tmpReal, tmpImag, RealOut, ImagOut);
    
    float wtemp = float (sin(0.5 * theta));
    
    float wpr = -2.0 * wtemp * wtemp;
    float wpi = float (sin(theta));
    float wr = 1.0 + wpr;
    float wi = wpi;
    
    int i3;
    
    float h1r, h1i, h2r, h2i, rt, it;
    //float total=0;
    
    for (i = 1; i < Half / 2; i++) {
        
        i3 = Half - i;
        
        h1r = 0.5 * (RealOut[i] + RealOut[i3]);
        h1i = 0.5 * (ImagOut[i] - ImagOut[i3]);
        h2r = 0.5 * (ImagOut[i] + ImagOut[i3]);
        h2i = -0.5 * (RealOut[i] - RealOut[i3]);
        
        rt = h1r + wr * h2r - wi * h2i; //printf("Realout%i = %f",i,rt);total+=fabs(rt);
        it = h1i + wr * h2i + wi * h2r; // printf("  Imageout%i = %f\n",i,it);
        
        Out[i] = rt * rt + it * it;
        
        rt = h1r - wr * h2r + wi * h2i;
        it = -h1i + wr * h2i + wi * h2r;
        
        Out[i3] = rt * rt + it * it;
        
        wr = (wtemp = wr) * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
    }
    //printf("total = %f\n",total);
    rt = (h1r = RealOut[0]) + ImagOut[0];
    it = h1r - ImagOut[0];
    Out[0] = rt * rt + it * it;
    
    rt = RealOut[Half / 2];
    it = ImagOut[Half / 2];
    Out[Half / 2] = rt * rt + it * it;
    
    delete[]tmpReal;
    delete[]tmpImag;
    delete[]RealOut;
    delete[]ImagOut;
}

/*
 * Windowing Functions
 */

int NumWindowFuncs()
{
    return 4;
}

char *WindowFuncName(int whichFunction)
{
    switch (whichFunction) {
        default:
        case 0:
            return "Rectangular";
        case 1:
            return "Bartlett";
        case 2:
            return "Hamming";
        case 3:
            return "Hanning";
    }
}

void WindowFunc(int whichFunction, int NumSamples, float *in)
{
    int i;
    
    if (whichFunction == 1) {
        // Bartlett (triangular) window
        for (i = 0; i < NumSamples / 2; i++) {
            in[i] *= (i / (float) (NumSamples / 2));
            in[i + (NumSamples / 2)] *=
            (1.0 - (i / (float) (NumSamples / 2)));
        }
    }
    
    if (whichFunction == 2) {
        // Hamming
        for (i = 0; i < NumSamples; i++)
            in[i] *= 0.54 - 0.46 * cos(2 * M_PI * i / (NumSamples - 1));
    }
    
    if (whichFunction == 3) {
        // Hanning
        for (i = 0; i < NumSamples; i++)
            in[i] *= 0.50 - 0.50 * cos(2 * M_PI * i / (NumSamples - 1));
    }
}

/* constructor */
fft::fft() {
    
}

/* destructor */
fft::~fft() {
    
    
}

/* Calculate the power spectrum */
void fft::powerSpectrum(int start, int half, float *data, int windowSize,float *magnitude,float *phase, float *power, float *avg_power) {
    int i;
    int windowFunc = 3;
    float total_power = 0.0f;
    
    /* processing variables*/
    float *in_real = new float[windowSize];
    float *in_img = new float[windowSize];
    float *out_real = new float[windowSize];
    float *out_img = new float[windowSize];
    
    for (i = 0; i < windowSize; i++) {
        in_real[i] = data[start + i];
    }
    
    WindowFunc(windowFunc, windowSize, in_real);
    RealFFT(windowSize, in_real, out_real, out_img);
    
    for (i = 0; i < half; i++) {
        /* compute power */
        power[i] = out_real[i]*out_real[i] + out_img[i]*out_img[i];
        total_power += power[i];
        /* compute magnitude and phase */
        magnitude[i] = 2.0*sqrt(power[i]);
        phase[i] = atan2(out_img[i],out_real[i]);
    }
    /* calculate average power */
    *(avg_power) = total_power / (float) half;
				
    delete[]in_real;
    delete[]in_img;
    delete[]out_real;
    delete[]out_img;
}

void fft::inversePowerSpectrum(int start, int half, int windowSize, float *finalOut,float *magnitude,float *phase) {
    int i;
    int windowFunc = 3;
    
    /* processing variables*/
    float *in_real = new float[windowSize];
    float *in_img = new float[windowSize];
    float *out_real = new float[windowSize];
    float *out_img = new float[windowSize];
    
    /* get real and imag part */
    for (i = 0; i < half; i++) {
        in_real[i] = magnitude[i]*cos(phase[i]);
        in_img[i]  = magnitude[i]*sin(phase[i]);
    }
    
    /* zero negative frequencies */
    for (i = half; i < windowSize; i++) {
        in_real[i] = 0.0;
        in_img[i] = 0.0;
    }
    
    FFT(windowSize, 1, in_real, in_img, out_real, out_img); // second parameter indicates inverse transform
    WindowFunc(windowFunc, windowSize, out_real);
				
    for (i = 0; i < windowSize; i++) {
        finalOut[start + i] += out_real[i];
    }
    
    delete[]in_real;
    delete[]in_img;
    delete[]out_real;
    delete[]out_img;
}

/*---------------後から追加-------------*/
void fft::setup(){
    band_bottom[0]=1,band_bottom[1]=3,band_bottom[2]=35,band_bottom[3]=60;
    band_top[0]=2,band_top[1]=10,band_top[2]=45,band_top[3]=70;
    for(int i=0;i<BAND_NUM;i++)map_min[i]=0;
    for(int i=0;i<BAND_NUM;i++)map_newMin[i]=0;
    map_max[0]=5,map_max[1]=1,map_max[2]=0.5,map_max[3]=0.5;
    map_newMax[0]=1,map_newMax[1]=1,map_newMax[2]=2,map_newMax[3]=2;
    rate=0.05;
    smoothRate = 0.7;
    for(int i=0;i<BAND_NUM;i++){
        lmh_length[i] = band_top[i] - band_bottom[i];
        pre_val[i]=1;
    }
}

void fft::update(float *magni,int i){
    for(int j=band_bottom[i];j<band_top[i];j++){
        if(magni[j]!=INFINITY)val[i]+=magni[j];
    }
    temp_val = val[i]/lmh_length[i];
    if(temp_val<map_min[i])bCut[i]=true;
    if(temp_val>vol_max[i]){
        vol_max[i]=temp_val;
        if(bAutoMaxGet)map_max[i]=vol_max[i];
    }
    
    //------マニュアルマップ関数
    float oldRate=temp_val/(map_max[i]-map_min[i]);
    float mapped = (map_newMax[i]-map_newMin[i])*oldRate+map_newMin[i];
    
    
    if(bSmooth){//平滑化あり
        pre_ave[i] = mapped;
        val[i] = smoothRate*pre_ave[i] + (1-smoothRate)*pre_val[i];
        pre_val[i]=val[i];
    }else{//平滑化なし
        val[i] = mapped;
    }
    
    if(bCut[i]){//一定値以下切り捨て
        bCut[i]=false;
        val[i]=0;
    }
    /*-------四捨五入----------*/
    float tmp = map_min[i]*100;
    tmp=round(tmp);
    map_min[i]=tmp/100;
    tmp = map_max[i]*100;
    tmp = round(tmp);
    map_max[i]=tmp/100;
}

void fft::changeBandRange(int key){
    if(key == 'q')band_bottom[0]++;
    else if(key=='a')band_bottom[0]--;
    else if(key=='w')band_top[0]++;
    else if(key=='s')band_top[0]--;
    else if(key=='e')band_bottom[1]++;
    else if(key=='d')band_bottom[1]--;
    else if(key=='r')band_top[1]++;
    else if(key=='f')band_top[1]--;
    else if(key=='t')band_bottom[2]++;
    else if(key=='g')band_bottom[2]--;
    else if(key=='y')band_top[2]++;
    else if(key=='h')band_top[2]--;
    else if(key=='u')band_bottom[3]++;
    else if(key=='j')band_bottom[3]--;
    else if(key=='i')band_top[3]++;
    else if(key=='k')band_top[3]--;
    for(int i=0;i<BAND_NUM;i++){
        lmh_length[i] = band_top[i] - band_bottom[i];
    }
}

void fft::changeParam(int key){
    if(key=='1')map_newMin[0]+=rate;
    else if(key=='q')map_newMin[0]-=rate;
    else if(key=='2')map_newMax[0]+=rate;
    else if(key=='w')map_newMax[0]-=rate;
    else if(key=='3')map_newMin[1]+=rate;
    else if(key=='e')map_newMin[1]-=rate;
    else if(key=='4')map_newMax[1]+=rate;
    else if(key=='r')map_newMax[1]-=rate;
    else if(key=='5')map_newMin[2]+=rate;
    else if(key=='t')map_newMin[2]-=rate;
    else if(key=='6')map_newMax[2]+=rate;
    else if(key=='y')map_newMax[2]-=rate;
    else if(key=='7')map_newMin[3]+=rate;
    else if(key=='u')map_newMin[3]-=rate;
    else if(key=='8')map_newMax[3]+=rate;
    else if(key=='i')map_newMax[3]-=rate;
    
    else if(key=='a')map_min[0]+=rate;
    else if(key=='z')map_min[0]-=rate;
    else if(key=='s')map_max[0]+=rate;
    else if(key=='x')map_max[0]-=rate;
    else if(key=='d')map_min[1]+=rate;
    else if(key=='c')map_min[1]-=rate;
    else if(key=='f')map_max[1]+=rate;
    else if(key=='v')map_max[1]-=rate;
    else if(key=='g')map_min[2]+=rate;
    else if(key=='b')map_min[2]-=rate;
    else if(key=='h')map_max[2]+=rate;
    else if(key=='n')map_max[2]-=rate;
    else if(key=='j')map_min[3]+=rate;
    else if(key=='m')map_min[3]-=rate;
    else if(key=='k')map_max[3]+=rate;
    else if(key==',')map_max[3]-=rate;
    
}
