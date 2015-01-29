/*
周波数重心を求めるプログラム
短時間フーリエ変換で求めた各時刻におけるスペクトルから周波数重心を求める
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include "wavelib.h"
#include "fftsg.h"
#define FFTLEN 4096
using namespace std;

void convertToDoubleFor8bit(unsigned char* org, double* conv);
void multHanningWindow(double* wav);
void executeFFT(double* wav, double* fft);
double computeCentroid(double* fft, double freqPerSample);

int main(void){
	struct RIFF riff;
	struct fmt fmt;
	struct data data;

	load_wave_data(&riff,&fmt,&data,"");

	if(fmt.BitsPerSample == 8 && fmt.Channel == 1){
		int length = data.size_of_sounds;
		int dataIndex = 0;
		unsigned char dividedWav[FFTLEN];
		double doubleWav[FFTLEN], centroid, currentTime;
		ofstream ofs("centroid.dat");

		double* fft = new double[2*FFTLEN];
		double freqPerSample = (double)fmt.SamplesPerSec / FFTLEN;

		while(length > dataIndex + FFTLEN){
			for(int i=0;i<FFTLEN;i++){
				dividedWav[i] = *(data.sounds + dataIndex + i);
			}
			convertToDoubleFor8bit(dividedWav,doubleWav);
			multHanningWindow(doubleWav);
			executeFFT(doubleWav, fft);
			centroid = computeCentroid(fft, freqPerSample);
			currentTime = (double)(dataIndex + FFTLEN/2) / fmt.SamplesPerSec;
			ofs << currentTime << "\t"
				<< centroid << endl;
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

double computeCentroid(double* fft, double freqPerSample){
	double power[FFTLEN/2];
	double powerSum = 0.0;
	double fpSum = 0.0;

	for(int i=0;i<FFTLEN/2;i++){
		power[i] = sqrt(fft[2*i]*fft[2*i] + fft[2*i+1]*fft[2*i+1]);
		powerSum += power[i];
		fpSum += (freqPerSample*i) * power[i];
	}
	return fpSum/powerSum;
}