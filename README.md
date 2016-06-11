//# TSP-
//中国34个城市的TSP问题。
// tsp.cpp : 定义控制台应用程序的入口点。
//
#include "stdafx.h"
#include <string.h>
#include <stdlib.h>
#include "math.h"
#include "time.h"
#define PI  3.1415926535898
#define R  6371.0;//km 地球半径 平均值，千米
#define CITY_NUM 34    
#define POPSIZE 100
#define MAXVALUE 10000   
#define N 10000
//unsigned seed = (unsigned)time(0);
int Hash[CITY_NUM + 1];
typedef struct {/*坐标结构*/
	char flag[10];/*标记城市*/
	float x;/*城市坐标*/
	float y;
}POSITION[CITY_NUM];
POSITION temp;
typedef char MAP[20][20];/*地图空间*/
double CityDistance[CITY_NUM][CITY_NUM];
MAP map;/*地图变量*/
typedef struct{
	int colony[POPSIZE][CITY_NUM + 1];
	double fitness[POPSIZE];
	double Distance[POPSIZE];//路径实际长度     针对不同的个体的路径总长度
	int BestRooting[CITY_NUM + 1];//最优城市路径序列
	double BestFitness;//最优路径适应值
	double BestValue;//最优路径长度
	int BestNum;
}TSP, *PTSP;
double ConvertDegreesToRadians(double degrees) {
	return degrees  * PI / 180;
}
double HaverSin(double theta) {
	double v = sin(theta / 2);
	return v*v;
}
double Distance(double wd1, double jd1, double wd2, double jd2) {// 根据经纬度坐标计算实际距离
	wd1 = ConvertDegreesToRadians(wd1);
	jd1 = ConvertDegreesToRadians(jd1);
	wd2 = ConvertDegreesToRadians(wd2);
	jd2 = ConvertDegreesToRadians(jd2);
	double wd = fabs(wd1 - wd2);
	double jd = fabs(jd1 - jd2);
	double h = HaverSin(wd) + cos(wd1)*cos(wd2)*HaverSin(jd);
	double distance = 2 * asin(sqrt(h)) * R;
	return distance;
}
void CalculatDist()
{			
	int i, j;
	int n = 1;
	//POSITION temp;
	for (i = 0;i < 34;i++)
	{
		//temp[i].flag ="空";
		CityDistance[i][i] = 0;
	}
	for (i = 0;i <20;i++) {
		for (j = 0;j < 20;j++)
		{
			map[i][j] = ' ';
		}
	}
	FILE *fp;
	if (NULL == (fp = fopen("shenghuijingweidu.txt", "r")))
	{
		printf("error\n");
		exit(1);
	}
	for (int i = 0;i < 34;i++)
	{
		fscanf(fp, "%s %f %f", &temp[i].flag, &temp[i].x, &temp[i].y);
		map[(int)(temp[i].y / 2) - 43][(int)(temp[i].x / 2) - 10] = '*';
		//printf("%s %3.2f %3.2f\n", temp[i].flag, temp[i].x, temp[i].y);//城市经纬度输出
	}
	fclose(fp);
	/*计算城市之间的距离词典*/
	for (i = 0;i < CITY_NUM;i++)
	{
		for (j = 0;j < CITY_NUM;j++)
		{
			if (i != j)
			{
				CityDistance[i][j] = Distance(temp[i].x, temp[i].y, temp[j].x, temp[j].y);
			}
		}
	}
	/*printf("城市地图\n\n");
	for (i = 0;i < 20;i++)
	{
		for (j = 0;j < 20;j++)
		{
			printf("%4c", map[i][j]);
		}
		printf("\n");
	}*/
	/*每一个城市的位置坐标*/
	/*for (i = 0;i < CITY_NUM;i++)
	{
		printf("%s的经纬度:lon=%3.2f,lat=%3.2f\n", temp[i].flag, temp[i].x, temp[i].y);
	}*/
	/*
	printf("\n城市距离词典:\n");
	for (i = 0;i < 34;i++)
	{
		printf("%s：  \n", temp[i].flag);
		for (j = 0;j < 34;j++)
		{
			printf("%lf   ", CityDistance[i][j]);
		}
		printf("\n\n");
	}*/
}
void copy(int a[], int b[])
{
	int i = 0;
	for (i = 0; i<CITY_NUM + 1; i++)
	{
		a[i] = b[i];
	}
}
bool check(TSP &city, int pop, int num, int k)
{//用来检查新生成的节点是否在当前群体中，0号节点是默认出发节点和终止节点
	int i;
	for (i = 0; i <= num; i++){
		if (k == city.colony[pop][i])
			return true;//新生成节点存在于已经生成的路径中
	}
	return false;//新生成节点没有存在于已经生成的路径中
}
void InitColony(TSP &city)
{
	int i, j, r;
	for (i = 0; i<POPSIZE; i++){
		city.colony[i][0] = 0;
		city.colony[i][CITY_NUM] = 0;
		city.BestValue = MAXVALUE;
		city.BestFitness = 0;//适应值越大越好
	}
	srand((unsigned)time(NULL));
	for (i = 0; i<POPSIZE; i++)
	{
		for (j = 1; j<CITY_NUM; j++)
		{
			r = rand() % (CITY_NUM - 1) + 1;//产生1～CITY_NUM-1之间的随机数 1-34 的随机数
			while (check(city, i, j, r))//一直产生城市节点
			{
				r = rand() % (CITY_NUM - 1) + 1;
			}
			city.colony[i][j] = r;
		}
	}
}
void CalFitness(TSP &city)//计算适应值
{
	int i, j;
	int start, end;
	int Best = 0;
	for (i = 0; i<POPSIZE; i++){//求适应值
		city.Distance[i] = 0;
		for (j = 1; j <= CITY_NUM; j++){
			start = city.colony[i][j - 1]; end = city.colony[i][j];
			city.Distance[i] = city.Distance[i] + CityDistance[start][end];
		}
		//	city.fitness[i]=pow(1.1,N/city.Distance[i]);
		city.fitness[i] = N / city.Distance[i];//个体适应值计算
		if (city.fitness[i]>city.fitness[Best])//求最大适应值也是个体的最优路径
			Best = i;
	}
	copy(city.BestRooting, city.colony[Best]);
	city.BestFitness = city.fitness[Best];
	city.BestValue = city.Distance[Best];
	city.BestNum = Best;//最优的个体值序号
}
void Select(TSP &city)
{//选择算子
	int TempColony[POPSIZE][CITY_NUM + 1];//[100][35]
	int i, s, t;
	double GaiLv[POPSIZE];
	int SelectP[POPSIZE + 1];
	double sum = 0;
	for (i = 0; i<POPSIZE; i++)
	{
		sum += city.fitness[i];
	}
	for (i = 0; i<POPSIZE; i++)
	{
		GaiLv[i] = city.fitness[i] / sum;
	}
	SelectP[0] = 0;
	for (i = 0; i<POPSIZE; i++)
	{
		SelectP[i + 1] = SelectP[i] + GaiLv[i] * RAND_MAX;//RAND_MAX是C中stdlib.h中宏定义的一个字符常量：#define RAND_MAX Ox7FFF其值最小为0, 最大为32767
		//	通常在产生随机小数时可以使用RAND_MAX。
	}
	memcpy(TempColony[0], city.colony[city.BestNum], sizeof(TempColony[0]));
	//	copy(TempColony[0],city.colony[city.BestNum]);
	for (t = 1; t<POPSIZE; t++)
	{
		s = rand() % RAND_MAX;//有待改进
		for (i = 1; i<POPSIZE; i++)
		{
			if (SelectP[i] >= s)
				break;
		}
		memcpy(TempColony[t], city.colony[i - 1], sizeof(TempColony[t]));
		//		copy(TempColony[t],city.colony[i-1]);
	}
	for (i = 0; i<POPSIZE; i++)
	{
		//		copy(city.colony[i],TempColony[i]);
		memcpy(city.colony[i], TempColony[i], sizeof(TempColony[i]));
	}
}
void Cross(TSP &city, double pc)
{//交叉概率是pc
	int i, j, t, l;
	int a, ca, cb;
	int Temp1[CITY_NUM + 1], Temp2[CITY_NUM + 1];
	for (i = 0; i<POPSIZE; i++)
	{
		double s = ((double)(rand() % RAND_MAX)) / RAND_MAX;
		if (s<pc)
		{
			ca = rand() % POPSIZE;
			cb = rand() % POPSIZE;
			while (cb != ca)cb = rand() % POPSIZE;
			if (ca == city.BestNum || cb == city.BestNum)
				continue;
			l = rand() % 17 + 1;
			a = rand() % (CITY_NUM - l) + 1;


			memset(Hash, 0, sizeof(Hash));
			Temp1[0] = Temp1[CITY_NUM] = 0;
			for (j = 1; j <= l; j++)
			{
				Temp1[j] = city.colony[cb][a + j - 1];
				Hash[Temp1[j]] = 1;
			}
			for (t = 1; t<CITY_NUM; t++)
			{
				if (Hash[city.colony[ca][t]] == 0)
				{
					Temp1[j++] = city.colony[ca][t];
					Hash[city.colony[ca][t]] = 1;
				}
			}
			memset(Hash, 0, sizeof(Hash));
			Temp2[0] = Temp2[CITY_NUM] = 0;
			for (j = 1; j <= l; j++)
			{
				Temp2[j] = city.colony[ca][a + j - 1];
				Hash[Temp2[j]] = 1;
			}
			for (t = 1; t<CITY_NUM; t++)
			{
				if (Hash[city.colony[cb][t]] == 0)//当Hash表内的某个城市为空时就是城市的值为0时就可以把cb的值赋给Temp2

				{
					Temp2[j++] = city.colony[cb][t];
					Hash[city.colony[cb][t]] = 1;
				}
			}
			//	copy(city.colony[ca],Temp1);
			//	copy(city.colony[cb],Temp2);
			memcpy(city.colony[ca], Temp1, sizeof(Temp1));
			memcpy(city.colony[cb], Temp2, sizeof(Temp2));
		}
	}
}
double GetFittness(int a[CITY_NUM + 1])
{
	int i, start, end;
	double Distance = 0;
	for (i = 0; i<CITY_NUM; i++)
	{
		start = a[i];   end = a[i + 1];
		Distance += CityDistance[start][end];
	}
	return N / Distance;
}
void Mutation(TSP &city, double pm)
{//变异概率是pm
	int i, k, m;
	int Temp[CITY_NUM + 1];
	for (k = 0; k<POPSIZE; k++)
	{
		double s = ((double)(rand() % RAND_MAX)) / RAND_MAX;
		i = rand() % POPSIZE;
		if (s<pm&&i != city.BestNum)
		{
			int a, b, t;
			a = (rand() % (CITY_NUM - 1)) + 1;
			b = (rand() % (CITY_NUM - 1)) + 1;
			copy(Temp, city.colony[i]);//i 是随机数，这样的话可以随机选择一个城市圈copy到Temp内
			if (a>b)
			{
				t = a;
				a = b;
				b = t;
			}
			for (m = a; m<(a + b) / 2; m++)
			{
				t = Temp[m];
				Temp[m] = Temp[a + b - m];
				Temp[a + b - m] = t;
			}
			if (GetFittness(Temp)<GetFittness(city.colony[i]))//是对其进行变异处理，把一个为i的城市圈进行内部交换城市元素，自身变异
			{
				a = (rand() % (CITY_NUM - 1)) + 1;
				b = (rand() % (CITY_NUM - 1)) + 1;
				//copy(Temp,city.colony[i]);
				memcpy(Temp, city.colony[i], sizeof(Temp));
				if (a>b)
				{
					t = a;
					a = b;
					b = t;
				}
				for (m = a; m<(a + b) / 2; m++)
				{
					t = Temp[m];
					Temp[m] = Temp[a + b - m];
					Temp[a + b - m] = t;
				}
			}
			memcpy(city.colony[i], Temp, sizeof(Temp));
		}
	}
}
void OutPut(TSP &city)
{
	int i;
	/*	printf("The population is:\n");
	for(i=0;i<POPSIZE;i++)
	{
	for(j=0;j<=CITY_NUM;j++)
	{
	printf("%5d",city.colony[i][j]);
	}
	printf("    %f\n",city.Distance[i]);
	}
	*/
	printf("最短的路径是:\n");
	for (i = 0; i <= CITY_NUM-1; i++)
	{
		printf("%s->", temp[(city.BestRooting[i])].flag);
	}
	printf("%s", temp[(city.BestRooting[CITY_NUM])].flag);
	printf("\n路径长度是:%f\n", (city.BestValue));
}
void Repeat(double pcross, double pmutation, int MaxEpoc) {
	TSP city;
	//srand(seed);
	CalculatDist();//求城市间两两之间的距离
				   //	city=(PTSP)malloc(sizeof(TSP));
	InitColony(city);//生成初始种群
	CalFitness(city);//计算适应值,考虑应该在这里面把最优选出来
	clock_t begin = clock();
	for (int i = 0; i<MaxEpoc; i++)
	{
		//	OutPut(city);
		Select(city);//选择(复制)
		Cross(city, pcross);//交叉     交叉率一般来说应该比较大，推荐使用80％-95％。
		Mutation(city, pmutation);//变异       变异率一般来说应该比较小，一般使用0.5％-1％最好。
		CalFitness(city);//计算适应值
	}
	clock_t end = clock();

	OutPut(city);//输出
	//printf("共用时间为：%lf秒\n", (double)(end - begin) / 1000);
	//system("pause");
}
int main()
{
	for (int i = 0;i < 100;i++) {
		printf("第%d次\n",i+1);
		Repeat((0.4+i/200), (0.095+i/20000), (3000+i*30));
	}
	return 0;
}
