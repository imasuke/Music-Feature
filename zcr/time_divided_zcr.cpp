/*
Zero-Crossing rate(ZCR)の時間変化を求めるプログラム
曲を25msで時分割し、分割した波形ごとにZCRを計算する
ここで計算するZCRはミリ秒単位の平均を取ったもの
*/
#include <iostream>
#include <fstream>
#include "wavelib.h"
#define DIVLEN 4096
using namespace std;

double compute_zcr_for8bit(unsigned char* wav, int length, int samples_per_sec);

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
		double zcr, currentTime;
		ofstream ofs("zcr.dat");

		while(length > data_index + divide_size){
			for(int i=0;i<divide_size;i++){
				divided_wav[i] = *(data.sounds + data_index + i);
			}
			zcr = compute_zcr_for8bit(divided_wav, divide_size, fmt.SamplesPerSec);
			currentTime = (double)(data_index + divide_size/2) / fmt.SamplesPerSec;
			ofs << currentTime << "\t" << zcr << endl;
			data_index += divide_size;
		}
	}
	else{
		cout << "このファイルには対応していません" << endl;
	}
	return 0;
}

double compute_zcr_for8bit(unsigned char* wav, int length, int samples_per_sec){
	unsigned char sample1, sample2;
	int zero_cross_count = 0;
	double zcr;

	for(int i=0;i<length-1;i++){
		sample1 = wav[i];
		sample2 = wav[i+1];
		if((sample1 > 128 && sample2 <= 128) || (sample1 < 128 && sample2 >= 128)){
			zero_cross_count++;
		}
	}
	zcr = (double)zero_cross_count / length * (samples_per_sec/1000);
	return zcr;
}