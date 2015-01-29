/*
Root Mean Square(RMS)を求めるプログラム
8bitモノラルwavファイルのみ対応
*/
#include<iostream>
#include<cmath>
#include"wavelib.h"
using namespace std;

int main(void){
	struct RIFF RIFF;
	struct fmt fmt;
	struct data data;
	load_wave_data(&RIFF,&fmt,&data,"");

	double square_mean = 0.0, wav;
	int length = data.size_of_sounds;
	double rms;
	if(fmt.Channel == 1 && fmt.BitsPerSample == 8){
		for(int i=0;i<length;i++){
			wav = (double)(*(data.sounds+i)-128);
			square_mean += (wav * wav) / length;
		}
		rms = sqrt(square_mean);
		cout << "rms = "<< rms << endl;
	}
	else{
		cout << "このファイルには対応していません" << endl;
	}
	return 0;
}