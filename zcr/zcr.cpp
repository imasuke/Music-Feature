/*
Zero-Crossing rate(ZCR)を求めるプログラム
1秒間のZCRを曲全体から平均して求める
8bitモノラルwavファイルのみ対応
*/
#include<iostream>
#include"wavelib.h"
using namespace std;

int main(void){
	struct RIFF RIFF;
	struct fmt fmt;
	struct data data;
	int zero_cross_count = 0;
	double zcr;

	load_wave_data(&RIFF,&fmt,&data,"");
	int length = data.size_of_sounds;
	unsigned char wav1,wav2;

	if(fmt.Channel == 1 && fmt.BitsPerSample == 8){
		for(int i=0;i<length-1;i++){
			wav1 = *(data.sounds+i);
			wav2 = *(data.sounds+i+1);
			if((wav1 > 128 && wav2 <= 128) || (wav1 < 128 && wav2 >= 128)){
				zero_cross_count++;
			}
		}
		zcr = (double)zero_cross_count / length * fmt.SamplesPerSec;
		cout << "zcr = " << zcr << endl;
	}
	else{
		cout << "このファイルには対応していません" << endl;
	}
	return 0;
}
