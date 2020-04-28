//#include <QCoreApplication>
#define _CRT_SECURE_NO_DEPRECATE；
#define _CRT_SECURE_NO_WARNINGS；
#pragma warning(disable:4996)；
#include <iostream>
#include <cstring>   
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "Gauss.h"
#include <fstream>    
#include <string>
#include <vector>
using namespace std;
#define PI          3.1415926535897932  /* pi */
#define D2R         (PI/180.0)          /* deg to rad */

/*
**	查找字符在字符串首次出现的位置
*/
int searchChar(char str[], char key) {
    //通过循环依次取得字符串每个字符，
    //循环结束条件 str[i] != '\0'
    for (int i = 0; str[i] != '\0'; i++) {
        if (str[i] == key) {
            return i;
        }
    }
    //如果没有结束，开始比较 key str[i] 
    return -1;
}
/*
**	高斯投影计算方法
*/
void Gauss(const double* BL, double* xy)
{
    //度转弧度
    double B = BL[1] * D2R; //纬度
    double L = BL[0] * D2R; //经度
    //椭球基本公式
    double a = 6378137;
    double f = 1 / 298.257222101;
    double L0 = 114;
    double b = a - a * f;
    double e2 = (a * a - b * b) / (a * a);
    double ep2 = (a * a - b * b) / (b * b);
    double M0 = a * (1 - e2);//子午圈赤道处的曲率半径
    //子午线弧长计算公式系数
    double Ac = 1 + 3 / 4.0 * e2 + 45 / 64.0 * pow(e2, 2) + 175 / 256.0 * pow(e2, 3) + 11025 / 16384.0 * pow(e2, 4) + 43659 / 65536.0 * pow(e2, 5);
    double Bc = 3 / 4.0 * e2 + 15 / 16.0 * pow(e2, 2) + 525 / 512.0 * pow(e2, 3) + 2205 / 2048.0 * pow(e2, 4) + 72765 / 65536.0 * pow(e2, 5);
    double Cc = 15 / 64.0 * pow(e2, 2) + 105 / 256.0 * pow(e2, 3) + 2205 / 4096.0 * pow(e2, 4) + 10395 / 16384.0 * pow(e2, 5);
    double Dc = 35 / 512.0 * pow(e2, 3) + 315 / 2048.0 * pow(e2, 4) + 31185 / 131072.0 * pow(e2, 5);
    double Ec = 315 / 16384.0 * pow(e2, 4) + 3465 / 65536.0 * pow(e2, 5);
    double Fc = 693 / 131072.0 * pow(e2, 5);

    double aerfa = Ac * M0;
    double beita = -1 / 2.0 * Bc * M0;
    double gama = 1 / 4.0 * Cc * M0;
    double delta = -1 / 6.0 * Dc * M0;
    double epsilon = 1 / 8.0 * Ec * M0;
    double zelta = -1 / 10.0 * Fc * M0;
    //大地坐标录入
    double W = sqrt(1 - e2 * sin(B) * sin(B));
    double N = a / W;
    double n2 = ep2 * cos(B) * cos(B);
    double t = tan(B);
    //子午线弧长
    double X = aerfa * B + beita * sin(2 * B) + gama * sin(4 * B) + delta * sin(6 * B) + epsilon * sin(8 * B) + zelta * sin(10 * B);
    double l = L - L0 * D2R; //经度差
    //辅助量计算
    double a0 = X;
    double a1 = N * cos(B);
    double a2 = 1 / 2.0 * N * pow(cos(B), 2) * t;
    double a3 = 1 / 6.0 * N * pow(cos(B), 3) * (1 - t * t + n2);
    double a4 = 1 / 24.0 * N * pow(cos(B), 4) * (5 - t * t + 9 * n2 + 4 * n2 * n2) * t;
    double a5 = 1 / 120.0 * N * pow(cos(B), 5) * (5 - 18 * t * t + pow(t, 4) + 14 * n2 - 58 * n2 * t * t);
    double a6 = 1 / 720.0 * N * pow(cos(B), 6) * (61 - 58 * t * t + pow(t, 4) + 270 * n2 - 330 * n2 * t * t) * t;
    //高斯正算
    xy[0] = a0 * pow(l, 0) + a2 * pow(l, 2) + a4 * pow(l, 4) + a6 * pow(l, 6); //x
    xy[1] = a1 * pow(l, 1) + a3 * pow(l, 3) + a5 * pow(l, 5); //y
    xy[1] += 500000;

}

int main()
{
    double* bl = new double[2], * xy = new double[2];
    int num = 0;
    FILE* fpRead = fopen("G:\\C-workspace\\Gauss\\省面.MIF", "r");     //需要修改
    FILE* fpWrite = fopen("G:\\C-workspace\\Gauss\\省面转换后.MIF", "w+");   //需要修改

    char StrLine[1024];   //每行最大读取的字符数
    if (fpRead == NULL)
    {
        printf("fpread open failed！");
        return -1;
    }
    if (fpWrite == NULL)
    {
        printf("fpWrite open failed！");
        return -1;
    }
    while (!std::feof(fpRead))
    {
        fgets(StrLine, 1024, fpRead);   //读取一行
        char* temp = new char[1024];
        if (strstr(StrLine, "Version") == NULL &&
            strstr(StrLine, "Charset") == NULL &&
            strstr(StrLine, "Delimiter") == NULL &&
            strstr(StrLine, "CoordSys") == NULL &&
            strstr(StrLine, "Columns") == NULL &&
            strstr(StrLine, "NAME") == NULL &&
            strstr(StrLine, "Data") == NULL &&
            strstr(StrLine, "Region") == NULL &&
            strstr(StrLine, "Pen") == NULL &&
            strstr(StrLine, "Brush") == NULL &&
            strstr(StrLine, "Center") == NULL &&
            strlen(StrLine) != 0)
        {
            if (strstr(StrLine, ".") == NULL)
            {
                strncpy(temp, StrLine, 10);
                bl[0] = atof(temp);
                xy[0] = bl[0];
                cout << xy[0] << endl;
            }
            else
            {
                int loc = searchChar(StrLine, ' ');
                strncpy(temp, StrLine, loc);
                bl[0] = atof(temp);
                strncpy(temp, StrLine+loc, strlen(StrLine)-loc);
                bl[1] = atof(temp);
                Gauss(bl, xy);
                cout << bl[0]<<","<<bl[1] << "," << xy[0] << "," << xy[1] << endl;
            }
            if ((int)bl[0] == 0)
            {
                continue;
            }
            if ((int)bl[0] == bl[0])
            {
                fprintf(fpWrite, "%13.f\n", bl[0]);
            }
            else
            {
                fprintf(fpWrite, "%13.7f %13.8f %13.7f %13.7f\n", bl[0], bl[1], xy[0], xy[1]);
            }
        }
    }
    bl = xy = NULL;
    std::fclose(fpRead);    //关闭文件
    std::fclose(fpWrite);   //关闭文件
    printf("END!\n");
    system("pause");
    return 0;
}
