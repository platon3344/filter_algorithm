#pragma once
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//% 程序功能  :采用JPDA数据关联算法实现两个个匀速运动目标的点迹与航迹的关联
//% 输入变量 :
//%           -target_position : 目标的初始位置
//%           -n : 采样次数
//%           -T : 采样间隔
//%           -MC_number : 仿真次数
//%           -c : 目标个数
//% 输出变量 :
//%           无
//% 参考文献 :
//%           黄玲, 数据挖掘及融合技术研究与应用, 西北工业大学硕士学位论文, 2004年
//% 声明      ：
//%           该代码为作者毕业设计内容，鉴于学术交流的角度，现在公开发布该代码
//%           该代码非本人原创，修改自网上另一位作者的JPDA代码
//%           该代码仅用于学术交流，请勿用于任何其它商业用途，请大家自觉遵守
//%           如果有人用该代码进行不合适的用途，该代码作者不承担任何责任
//%           请遵守作者的劳动成果，转载请标明
//% 作者邮箱  ：
//%           wangzexun@gmail.com
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#include <Dense>
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>

namespace JPDA{
using namespace std;
using namespace Eigen;

const double gc_dPi = 3.1415926;

class JPDA
{
public:
	JPDA(){
		init();
	}
	~JPDA(){}

public:
	void test();

	void init();

	void setEstimate(MatrixXd data){
		for (int i = 0; i < data.cols(); i++){
			xEstimate.at(i) = data.col(i);
		}
	}
	void getEstimate(vector<MatrixXd>& data){
		data = xEstimate;
	}
	void update(MatrixXd measureData);
	static void readTxtStreamData(string filename, vector<MatrixXd>& data);
	static void writeTxtStreamData(MatrixXd data, string filename);

private:
	//初始化
	int initFlag;
	//double 目标个数，在真实的环境中，应该是航迹的个数
	double m_iTargetNum;
	//周期
	double m_dPeriod;
	//检测概率
	double pd;
	//正确量测落入跟踪门的概率
	double pg;
	//门限
	double gSigma;
	//每个单位面积（km^2）内产生lambda个杂波
	double lambda;
	double gamma;
	//测量值误差
	MatrixXd targetDelta;
	//协方差矩阵
	vector<MatrixXd> p;
	//状态转移矩阵
	MatrixXd phi;
	//观测矩阵
	MatrixXd h;
	//观测协方差矩阵R
	MatrixXd R;
	//过程噪声协方差矩阵Q
	MatrixXd Q;
	//过程噪声矩阵G
	MatrixXd G;
	//目标各个时刻的滤波值
	vector<MatrixXd> xFilter;

	//多个目标，椭圆跟踪门面积
	vector<double> ellipseVol;
	//多个目标，中间变量
	vector<MatrixXd> S;
	//多个目标，预测值(x,vx,y,vy)
	vector<MatrixXd> xPredict;
	//多个目标，估计值坐标(x,vx,y,vy)
	vector<MatrixXd> xEstimate;
	//多个目标，协方差矩阵
	vector<MatrixXd> pPredict;
	//多个目标，预测坐标值(x,y)
	vector<MatrixXd> zPredict;
	//
};
}

