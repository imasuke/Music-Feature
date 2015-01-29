/*
周波数重心を求めるプログラム
入力したwavファイルの中心部分を抽出し、その周波数重心を求める
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include "wavelib.h"
#include "fftsg.h"
#define FFTLEN 4096
using namespace std;

void convert_to_double_for_16bit(unsigned char* org, double* conv);
void mult_hanning_window(double* wav);
void execute_FFT(double* wav, double* fft);
double compute_centroid(double* fft, double freqPerSample);

int main(int argc, char *argv[]){
	struct RIFF riff;
	struct fmt fmt;
	struct data data;

	//引数のwavファイルを読み込む
	load_wave_data(&riff,&fmt,&data,argv[1]);

	//16bitモノラルならば処理を続ける
	if(fmt.BitsPerSample == 16 && fmt.Channel == 1){
		double centroid, freq_per_sample;

		//配列の確保
		unsigned char *divided_wav = new unsigned char[FFTLEN*2];
		double *double_wav = new double[FFTLEN];
		double *fft = new double[2*FFTLEN];

		//中心のサンプルを抽出する
		int center = data.size_of_sounds/2;
		int st = center - FFTLEN;
		for(int i=0;i<FFTLEN*2;i++){
			divided_wav[i] = data.sounds[st+i];
		}
		//unsigned charをdoubleに変換
		convert_to_double_for_16bit(divided_wav,double_wav);
		//ハニング窓をかける
		mult_hanning_window(double_wav);
		//FFTを実行
		execute_FFT(double_wav,fft);
		//1サンプル当たりの周波数幅を求める
		freq_per_sample = (double)fmt.SamplesPerSec / FFTLEN;
		//周波数重心を求める
		centroid = compute_centroid(fft, freq_per_sample);
		//標準出力へ周波数重心を出力
		cout << centroid << endl;

		//メモリの解放
		delete[] divided_wav;
		delete[] double_wav;
		delete[] fft;
	}
	else{
		cerr << "このファイルには対応していません" << endl;
	}
	return 0;
}

void convert_to_double_for_16bit(unsigned char* org, double* conv){
	short s;
	for(int i=0;i<FFTLEN;i++){
		s = (org[2*i+1] << 8) + org[2*i];
		conv[i] = (double)s;
	}
}

void mult_hanning_window(double* wav){
	for(int i=0;i<FFTLEN;i++){
		wav[i] *=  0.5 - 0.5 * cos(2*M_PI*((double)i/FFTLEN));
	}
}

void execute_FFT(double* wav, double* fft){
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

double compute_centroid(double* fft, double freq_per_sample){
	double power[FFTLEN/2];
	double power_sum = 0.0;
	double fp_sum = 0.0;

	for(int i=0;i<FFTLEN/2;i++){
		power[i] = sqrt(fft[2*i]*fft[2*i] + fft[2*i+1]*fft[2*i+1]);
		power_sum += power[i];
		fp_sum += (freq_per_sample*i) * power[i];
	}
	return fp_sum/power_sum;
}