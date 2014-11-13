#include "stdio.h"
#include "stdlib.h"
#define N 1024
#define Pi 3.1415926
#include "math.h"
//DFT算法
void DFT(float *Input, int Amount, float *outputRe, float *outputIm){
/*Input 输入数组，Amount 数组长度 outputRe 实数输出部分 outputIm虚数部分 */
	int n, k;//循环控制变量
	float Cos, Sin;//临时变量
/*
	X[k]=x[n]*e^-i*2pi/N*k
	e^-i*2pi/N*k 欧拉变换[e^-x=cos(x)-isin(x)]
*/
	for (k = 0; k < Amount ;k++)
	{

		for (n = 0; n < Amount ; n++)
		{
			Cos = cos(2 * Pi*n*k / Amount);
			Sin =-sin(2 * Pi*n*k / Amount);
			outputRe[k] += Input[n] * Cos;
			outputIm[k] += Input[n] * Sin;
		}
		
	}

}
//递归实现FFT
void  FFT(float *Input, int Amount, float *outputRe, float *outputIm){
	float Wave1[N/2];
	float Wave2[N/2];//Wave1/2数组用于分别存放输入数组的奇数部分，偶数部分
	float C1;
	float S1;//cos与sin结果的临时变量
	int HF = Amount / 2,i;
	float CLW1_Re[N] = {0}, CLW1_Im[N] = {0}, CLW2_Re[N] = { 0 }, CLW2_Im[N] = { 0 };//W1,W2临时结果数组
	
	for ( i= 0; i < HF; i++)
	{
		Wave1[i] = Input[i * 2], Wave2[i] = Input[1 + i * 2];//分解输入数组
	}
	//当数组长度小于16时调用DFT计算结果
	if (Amount<16)
	{
		DFT(Wave1, HF, CLW1_Re, CLW1_Im);
		DFT(Wave2, HF, CLW2_Re, CLW2_Im);
	}
	else//数组长度大于16进入FFT递归分解
	{
		FFT(Wave1, HF, CLW1_Re, CLW1_Im);
		FFT(Wave2, HF, CLW2_Re, CLW2_Im);
	}
	//DFT/FFT的计算结果实数部分前后对称
	//虚数部分正负对称
	for ( i = 0; i < HF; i++)
		{
			CLW1_Re[i + HF] = CLW1_Re[HF-i];
			CLW2_Re[i + HF] = CLW2_Re[HF - i];
			CLW1_Im[i + HF] = -CLW1_Im[HF - i];
			CLW2_Im[i + HF] = -CLW2_Im[HF - i];
		}
	//FFT[x][k]=DFT[x1][k]*DET[x2][k]*(cos(2*pi/n*k)-isin(2*pi/n*k))
	for ( i = 0; i < Amount; i++)
	{
		C1=cos(2 * Pi / Amount*i) ;
		S1=-sin(2 * Pi / Amount*i) ;
		outputRe[i] = CLW1_Re[i] +C1 * CLW2_Re[i] -S1 *CLW2_Im[i];
		outputIm[i] = CLW1_Im[i] +S1 *CLW2_Re[i] +C1*CLW2_Re[i];
	}

}

/*DFT/FFT的计算结果实数部分前后对称粗略推导
X[3] += x[n] * [cos(2Pi/N * n * 3) - isin(2Pi/N * n * 3)];

X[N - 3] += x[n] * [cos(2Pi/N * n * (N - 3)) - isin(2Pi/N * n * (N - 3))]
=>X[N - 3] += x[n] * [cos(2Pi * n - 3 * 2Pi/N * n) - isin(2Pi * n - 3 * 2Pi/N * n)]
=>X[N - 3] += x[n] * [cos(- 3 * 2Pi/N * n) - isin(- 3 * 2Pi/N * n)]
由于cos(-x) = cos(x)，sin(-x) = - sin(x)
=>X[N - 3] += x[n] * [cos(3 * 2Pi/N * n) + isin(3 * 2Pi/N * n)]
所以
X[N - 3]的实数部分 = X[3]的实数部分
X[N - 3]的虚数部分 = - X[3]的虚数部分
*/

//碟型FFT FFT(Wave)[k]=DFT(x1)[k]+WKN*DFT(x2)[k];
int  RxFFT(int Power, float *outputRe, float *outputIm){
	float TPQA,TPQATP;//存放WKN的计算值
	float Re,Im,ReFactor,ImFactor;//各项乘法加法临时存放变量
	int i,j,k,antother,current,currentBase;
	//i，k三层循环的控制变量 current,antother碟型变换对应两个变量，currentBase计算交换变量的中介值
	int Amount,Layer,pmul;
	//Amount数组长度，Layer层数，j交换的变量数,pmul每一个小型fft的个数

	if (Power<1||Power>31)//防止溢出
	{
		return 0;
	}
	Amount=pow(2,Power);//获取数组长度
	TPQA=Pi*2/Amount;
	for (Layer = 1; Layer <= Power; ++Layer)
	{
		j=pow(2,Layer-1);//2^layer-1
		pmul=pow(2,Power-Layer);//计算小型fft的个数
		for (i = 0; i < pmul; ++i)
		{
			currentBase=i*2*j;
			TPQATP=TPQA*pmul;
			//第一层的碟型计算
			{
				current=currentBase;
				antother=current+j;
				//碟型加法
				Re=outputRe[current]+outputRe[antother];
				Im=outputIm[current]+outputIm[antother];
				//因为是一层WNK分别等于=1跟0，所以直接交换
				outputRe[antother]=outputRe[current]-outputRe[antother];
				outputIm[antother]=outputIm[current]-outputIm[antother];
				//
				outputRe[current]=Re;
				outputIm[current]=Im;
			}

			//每一层碟型运算循环
			for ( k = 1; k < j; ++k)
			{
				current=currentBase+k;
				antother=current+j;
				//计算WNK
				ReFactor=cos(TPQATP*k);
				ImFactor=-sin(TPQATP*k);
				//乘与WNK
				Re=outputRe[antother];
				outputRe[antother]=outputRe[antother]*ReFactor - outputIm[antother]*ImFactor;
				outputIm[antother]=outputIm[antother]*ReFactor+Re*ImFactor;
				//碟型加法
				Re=outputRe[current];
				Im=outputIm[current];
				outputRe[current]+=outputRe[antother];
				outputIm[current]+=outputIm[antother];
				//
				outputRe[antother]=Re-outputRe[antother];
				outputIm[antother]=Im-outputIm[antother];
			}
		}
	}
	return 1;
}
//乱序排序
/*
例子
1=001=>100=4
2=010=>010=2
3=011=>110=6
.....

*/
int REversrArrang_Auto(float NewArray[],float Array[],int power) {
	int Amount, i,j;
	int tmp=0x00000000;
	if(power<1||power>31)
		return 0;
	Amount=pow(2,power);
	NewArray[0]=Array[0];
	//从0到FFT大小-1的循环
	for (i = 0; i < Amount-1; ++i)
	{		
		j=power-1;//检查起始数位
//逆向二进制加法
		while((tmp &(1<<j))!=0)//直到地j位为零
		{
			tmp&=~(1<<j);//把第j位归零
			j--;//跳到前一位
		}
		tmp|= (1<<j);//把第j位算一
		NewArray[i+1]=Array[tmp];//按新顺序输出
	}
	return 1;
}

//碟型原型
int Radix2_Dynamic(float Real[], float Imag[], int Power)
{
	int Layer;
	int i, j, k, pmul, current, another, currentBase;
	int Amount;
	float ReFactor, ImFactor;
	float Re, Im;
	float TPOA, TPOATP; //TPOA = TwoPiOverAmount, TPOATP = TwoPiOverAmountTimesPmul
	if(Power > 31 || Power < 1)
		return 0;
	Amount = pow(2, Power);
	TPOA = 2 * 3.1415926535 / Amount;

	//For each layer in fft.
    for(Layer = 1; Layer <= Power; Layer ++)
	{
		j = pow(2, Layer - 1);
		pmul = pow(2, Power - Layer); //W_Amount ^ (0 to p interval j)
		
        for(i = 0; i < pmul; i ++)
		{
			currentBase = i * 2 * j;
			TPOATP = TPOA * pmul;
			{
				//k = 0
				current = currentBase;
				another = current + j;
				
				Re = Real[current] + Real[another];
				Im = Imag[current] + Imag[another];
				Real[another] = Real[current] - Real[another];
				Imag[another] = Imag[current] - Imag[another];
				
				Real[current] = Re;
				Imag[current] = Im;
			}
			{
                for(k = 1; k < j; k ++)
				{
					current = currentBase + k;
					another = current + j;

					
						ReFactor = cos(TPOATP * k);
						ImFactor = -sin(TPOATP * k);

						Re = Real[another];
						Real[another] = Real[another] * ReFactor - Imag[another] * ImFactor;
						Imag[another] = Imag[another] * ReFactor + Re * ImFactor;

						Re = Real[current];
						Im = Imag[current];
						Real[current] += Real[another];
						Imag[current] += Imag[another];

						Real[another] = Re - Real[another];
						Imag[another] = Im - Imag[another];
					
				}
			}
		}
	}
	return 1;
}
//主测试函数
int main(int argc, char const *argv[]){
	float inp[N] = { 0};
	float Re[N] = { 0 };
	float Im[N] = { 0 };
	float Re1[N] = { 0 };
	float Im1[N] = { 0 };
	int i;
	for ( i = 0; i < N; i++)
	{ 
		inp[i] = i;
	}
	
         /*	FFT(inp, N, Re1, Im1);

	for ( i = 0; i < N; i++)
	{
		printf("%0.2f+%0.2f\n", Re1[i], Im1[i]);
	}*/
	REversrArrang_Auto(Re,inp,10);
	//RxFFT(10,Re,Im);
	Radix2_Dynamic(Re,Im,10);
	for ( i = 0; i < N; i++)
	{
		printf("%f+%0.2f\n", Re[i], Im[i]);
	}
	system("pause");
	return 0;
}
