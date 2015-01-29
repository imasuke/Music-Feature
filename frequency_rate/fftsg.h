#ifndef FFTSG_H
#define FFTSG_H

/*複素DFT
n    : FFT点列*2の数
isgn : -1で離散フーリエ変換,1で離散フーリエ逆変換
a    : 入出力配列. 長さはFFT点列の2倍. 偶数番目を実部,奇数番目を虚部とする.
ip   : 作業用配列. 長さは2+sqrt(n)以上必要. 最初にcdft()を呼ぶ時にはip[0]=0を代入しておく.
w    : 三角関数用配列. 長さはn/2とする.
*/
void cdft(int n, int isgn, double *a, int *ip, double *w);

void rdft(int n, int isgn, double *a, int *ip, double *w);
void ddct(int n, int isgn, double *a, int *ip, double *w);
void ddst(int n, int isgn, double *a, int *ip, double *w);
void dfct(int n, double *a, double *t, int *ip, double *w);
void dfst(int n, double *a, double *t, int *ip, double *w);

#endif
