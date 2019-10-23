#pragma once
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//% ������  :����JPDA���ݹ����㷨ʵ�������������˶�Ŀ��ĵ㼣�뺽���Ĺ���
//% ������� :
//%           -target_position : Ŀ��ĳ�ʼλ��
//%           -n : ��������
//%           -T : �������
//%           -MC_number : �������
//%           -c : Ŀ�����
//% ������� :
//%           ��
//% �ο����� :
//%           ����, �����ھ��ںϼ����о���Ӧ��, ������ҵ��ѧ˶ʿѧλ����, 2004��
//% ����      ��
//%           �ô���Ϊ���߱�ҵ������ݣ�����ѧ�������ĽǶȣ����ڹ��������ô���
//%           �ô���Ǳ���ԭ�����޸���������һλ���ߵ�JPDA����
//%           �ô��������ѧ�����������������κ�������ҵ��;�������Ծ�����
//%           ��������øô�����в����ʵ���;���ô������߲��е��κ�����
//%           ���������ߵ��Ͷ��ɹ���ת�������
//% ��������  ��
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
	//��ʼ��
	int initFlag;
	//double Ŀ�����������ʵ�Ļ����У�Ӧ���Ǻ����ĸ���
	double m_iTargetNum;
	//����
	double m_dPeriod;
	//������
	double pd;
	//��ȷ������������ŵĸ���
	double pg;
	//����
	double gSigma;
	//ÿ����λ�����km^2���ڲ���lambda���Ӳ�
	double lambda;
	double gamma;
	//����ֵ���
	MatrixXd targetDelta;
	//Э�������
	vector<MatrixXd> p;
	//״̬ת�ƾ���
	MatrixXd phi;
	//�۲����
	MatrixXd h;
	//�۲�Э�������R
	MatrixXd R;
	//��������Э�������Q
	MatrixXd Q;
	//������������G
	MatrixXd G;
	//Ŀ�����ʱ�̵��˲�ֵ
	vector<MatrixXd> xFilter;

	//���Ŀ�꣬��Բ���������
	vector<double> ellipseVol;
	//���Ŀ�꣬�м����
	vector<MatrixXd> S;
	//���Ŀ�꣬Ԥ��ֵ(x,vx,y,vy)
	vector<MatrixXd> xPredict;
	//���Ŀ�꣬����ֵ����(x,vx,y,vy)
	vector<MatrixXd> xEstimate;
	//���Ŀ�꣬Э�������
	vector<MatrixXd> pPredict;
	//���Ŀ�꣬Ԥ������ֵ(x,y)
	vector<MatrixXd> zPredict;
	//
};
}

