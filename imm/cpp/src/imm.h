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
	
	//������immģ��
	int m_iIMMState = 6;
	int m_iModelNum = 3;
	vector<MatrixXd> imm_xn_est;
	vector<MatrixXd> imm_pn_est;
	MatrixXd imm_x_est_exd, imm_p_est_exd;	//����������λ��x_est_exd
	//����ģ�͸���ת������
	MatrixXd immU;
	//����ģ��ת��������Ʒ����ת�ƾ���
	MatrixXd immP;
	//����������ͬ��ģ�Ͳ�����ģ��1���ǻ���ģ�ͣ�ģ��2��3Ϊ����ģ�ͣ�Qֵ��ͬ
	vector<MatrixXd> immPhi;
	//��������Э�������g
	vector<MatrixXd> immg;
	//����������������g, Q = immg * immq * immg;
	vector<MatrixXd> immq;
	//��������
	MatrixXd immH;
	//������������ R
	MatrixXd immR;
	//��ϸ��ʾ���
	MatrixXd immMu;
	//��һ������
	MatrixXd immCmean;

};
}
