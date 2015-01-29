/*
Root Mean Square(RMS)の時間変化を求めるプログラム
曲を25msで時分割し、分割した波形ごとにRMSを計算する
*/
#include <iostream>
#include <cmath>
#include <fstream>
#include "wavelib.h"
#define DIVLEN 4096 
using namespace std;

double compute_rms_for8bit(unsigned char* wav, int length);

int main(void){
	struct RIFF RIFF;
	struct fmt fmt;
	struct data data;

	load_wave_data(&RIFF,&fmt,&data,"");

	if(fmt.BitsPerSample == 8 && fmt.Channel == 1){
		int length = data.size_of_sounds;
		int divide_size = DIVLEN;
		unsigned char divided_wav[divide_size];
		int data_index = 0;
		double rms, currentTime;
		ofstream ofs("rms.dat");

		while(length > data_index + divide_size){
			for(int i=0;i<divide_size;i++){
				divided_wav[i] = *(data.sounds + data_index + i);
			}
			rms = compute_rms_for8bit(divided_wav, divide_size);
			currentTime = (double)(data_index + divide_size/2) / fmt.SamplesPerSec;
			ofs << currentTime << "\t" << rms << endl;
			data_index += divide_size;
		}
	}
	else{
		cout << "このファイルには対応していません" << endl;
	}
	return 0;
}

double compute_rms_for8bit(unsigned char* wav, int length){
	double square_mean = 0.0, sample, rms;

	for(int i=0;i<length;i++){
		sample = (double)wav[i] - 128;
		square_mean +=  sample * sample;
	}
	rms = sqrt(square_mean/length);
	return rms;
}