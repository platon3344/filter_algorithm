#pragma once
#include <Dense>
#include <iostream>
#include <fstream>
#include <vector>

namespace IMM{
using namespace std;
using namespace Eigen;
const double gc_dPi = 3.1415926;

class imm
{
public:
	imm(){
		initKF();
		initIMM();
	}
	~imm(){}

public:
	void initKF();
	void initIMM();
	void test();
	void setEstimateData(double x, double y){
		x_est.resize(4, 1);
		x_est << x, 0,y,0;
	}
	void setEstimateData(double x, double vx, double y, double vy){
		x_est.resize(4, 1);
		x_est << x, vx, y, vy;
	}
	void setMeasureData(double x, double y){
		z.resize(2, 1);
		z << x, y;
	}
	void updateIMM();
	void updateKF();

	void getEstData(MatrixXd& data){
		data << imm_x_est_exd(0, 0) , imm_x_est_exd(0, 2);
	}

	static void readTestData(string filename, MatrixXd& data);
	static void writeResultData(MatrixXd data, string filename);

private:
	MatrixXd phi, h, g, R,z;
	MatrixXd x_est, p_est;
	double m_iPeriod = 2;
	double RDelta;
	int m_iKalmanState = 4;
	
	//以下是imm模型
	int m_iIMMState = 6;
	int m_iModelNum = 3;
	vector<MatrixXd> imm_xn_est;
	vector<MatrixXd> imm_pn_est;
	MatrixXd imm_x_est_exd, imm_p_est_exd;	//输出结果保存位置x_est_exd
	//三个模型概率转换矩阵
	MatrixXd immU;
	//控制模型转换的马尔科夫概率转移矩阵
	MatrixXd immP;
	//采用三个不同的模型参数，模型1个非机动模型，模型2、3为机动模型（Q值不同
	vector<MatrixXd> immPhi;
	//过程噪声协方差矩阵g
	vector<MatrixXd> immg;
	//过程噪声噪声矩阵g, Q = immg * immq * immg;
	vector<MatrixXd> immq;
	//测量矩阵
	MatrixXd immH;
	//测量噪声矩阵 R
	MatrixXd immR;
	//混合概率矩阵
	MatrixXd immMu;
	//归一化常数
	MatrixXd immCmean;

};
}
