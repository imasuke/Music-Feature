/*
短時間フーリエ変換を行ない、周波数分布に占める各周波数帯域の割合を求めるプログラム
割合を求める際には4つの周波数帯域に分ける
low... 0 <= f< 4000
middle... 4000 <= f < 8000
high... 8000 <= f < 12000
ultrahigh... 12000 <= f
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include "wavelib.h"
#include "fftsg.h"
#define FFTLEN 4096
#define LOW 4000 
#define MIDDLE 8000 
#define HIGH 12000
#define BANDNUM 4
using namespace std;

void convertToDoubleFor8bit(unsigned char* org, double* conv);
void multHanningWindow(double* wav);
void executeFFT(double* wav, double* fft);
void computeFreqRate(double* fft, double* freqRate, double freqPerSample);

int main(void){
	struct RIFF riff;
	struct fmt fmt;
	struct data data;

	load_wave_data(&riff,&fmt,&data,"");

	if(fmt.BitsPerSample == 8 && fmt.Channel == 1){
		int length = data.size_of_sounds;
		int dataIndex = 0;
		unsigned char dividedWav[FFTLEN];
		double doubleWav[FFTLEN], freqRate[BANDNUM], currentTime;
		ofstream ofs("freq_rate.dat");

		double* fft = new double[2*FFTLEN];

		while(length > dataIndex + FFTLEN){
			for(int i=0;i<FFTLEN;i++){
				dividedWav[i] = *(data.sounds + dataIndex + i);
			}
			convertToDoubleFor8bit(dividedWav,doubleWav);
			multHanningWindow(doubleWav);
			executeFFT(doubleWav, fft);
			computeFreqRate(fft, freqRate, (double)fmt.SamplesPerSec/FFTLEN);
			currentTime = (double)(dataIndex + FFTLEN/2) / fmt.SamplesPerSec;
			ofs << currentTime << "\t"
				<< freqRate[0] << "\t" 
				<< freqRate[1] << "\t"
				<< freqRate[2] << "\t"
				<< freqRate[3] << endl;
			dataIndex += FFTLEN;
		}
		delete[] fft;
	}
	else{
		cout << "このファイルには対応していません" << endl;
	}
	return 0;
}

void convertToDoubleFor8bit(unsigned char* org, double* conv){
	for(int i=0;i<FFTLEN;i++){
		conv[i] = (double)org[i] - 128;
	}
}

void multHanningWindow(double* wav){
	for(int i=0;i<FFTLEN;i++){
		wav[i] *=  0.5 - 0.5 * cos(2*M_PI*((double)i/FFTLEN));
	}
}

void executeFFT(double* wav, double* fft){
	int* ip = new int[2+(int)sqrt(2*FFTLEN)];
	double* w = new double[FFTLEN];

	for(int i=0;i<FFTLEN;i++){
		fft[2*i] = wav[i];
		fft[2*i+1] = 0.0;
	}
	ip[0] = 0;
	cdft(FFTLEN*2, -1, fft, ip, w);
	delete[] ip;
	delete[] w;
}

void computeFreqRate(double* fft, double* freqRate, double freqPerSample){
	double power[FFTLEN/2], bandSum[BANDNUM];
	double powerSum = 0.0, currentFreq;

	for(int i=0;i<BANDNUM;i++){
		bandSum[i] = 0.0;
	}
	for(int i=0;i<FFTLEN/2;i++){
		/*振幅スペクトルの計算*/
		power[i] = sqrt(fft[2*i]*fft[2*i] + fft[2*i+1]*fft[2*i+1]);
		if(power[i] < 1){
			power[i] = 1;
		}
		power[i] = 10*log10(power[i]);

		powerSum += power[i];
		currentFreq = i * freqPerSample;
		if(currentFreq < LOW){
			bandSum[0] += power[i];
		}
		else if(currentFreq < MIDDLE){
			bandSum[1] += power[i];
		}
		else if(currentFreq < HIGH){
			bandSum[2] += power[i];
		}
		else{
			bandSum[3] += power[i];
		}
	}
	for(int i=0;i<BANDNUM;i++){
		freqRate[i] = bandSum[i] / powerSum;
	}
}
