#include "imm.h"

namespace IMM{
void imm::initKF()
{
	//1 �����Ǳ�׼kalman�˲�ģ�ͣ�����4άģ�ͣ�x,vx,y,vy;
	//�˶�ģ��
	RDelta = 100.0;//�����������,��ֵ��matlab�������ݸ���
	phi.resize(4, 4);
	phi << 1, m_iPeriod, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, m_iPeriod,
		0, 0, 0, 1;

	//��������
	h.resize(2, 4);
	h << 1, 0, 0, 0,
		0, 0, 1, 0;

	g.resize(4, 2);
	g << m_iPeriod / 2, 0,
		1, 0,
		0, m_iPeriod,
		0, 1;
	//������������
	R.resize(2, 2);
	R << RDelta*RDelta, 0,
		0, RDelta*RDelta;
	z.resize(2, 1);//����ֵ

	x_est.resize(4, 1);
	p_est.resize(4, 4);

	p_est << RDelta*RDelta, RDelta*RDelta / m_iPeriod, 0, 0,
		RDelta*RDelta / m_iPeriod, 2 * RDelta*RDelta / (m_iPeriod*m_iPeriod), 0, 0,
		0, 0, RDelta*RDelta, RDelta*RDelta / m_iPeriod,
		0, 0, RDelta*RDelta / m_iPeriod, 2 * RDelta*RDelta / (m_iPeriod*m_iPeriod);
}
void imm::initIMM()
{
	//2 ����Ϊimm�˲�ģ��
	for (int i = 0; i < m_iModelNum; i++){
		MatrixXd tmp;
		tmp.resize(m_iIMMState, m_iIMMState);
		imm_pn_est.push_back(tmp);//(6,6,3)

		MatrixXd tmpx;
		tmpx.resize(m_iIMMState, 1);
		imm_xn_est.push_back(tmpx);//(6,1,3)
	}
	immU.resize(1, m_iModelNum);
	immU << 0.8, 0.1, 0.1;

	imm_x_est_exd.resize(m_iIMMState, 1);
	imm_x_est_exd.setZero();
	imm_p_est_exd.resize(m_iIMMState, m_iIMMState);
	imm_p_est_exd.setZero();

	immP.resize(m_iModelNum, m_iModelNum);
	//ÿһ�д���ǰ��ģ�͵�ת�Ƹ���
	immP << 0.95, 0.025, 0.025,
		0.025, 0.95, 0.025,
		0.025, 0.025, 0.95;
	MatrixXd tmp;
	tmp.resize(m_iIMMState, m_iIMMState);
	tmp.setZero();
	tmp << 1, m_iPeriod, 0, 0, 0, 0,
		0, 1, 0, 0, 0, 0,
		0, 0, 1, m_iPeriod, 0, 0,
		0, 0, 0, 1, 0, 0,
		0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0;
	immPhi.push_back(tmp);
	tmp << 1, m_iPeriod, 0, 0, m_iPeriod*m_iPeriod / 2, 0,
		0, 1, 0, 0, m_iPeriod, 0,
		0, 0, 1, m_iPeriod, 0, m_iPeriod*m_iPeriod / 2,
		0, 0, 0, 1, 0, m_iPeriod,
		0, 0, 0, 0, 1, 0,
		0, 0, 0, 0, 0, 1;
	immPhi.push_back(tmp);
	immPhi.push_back(tmp);

	MatrixXd gTmp;
	gTmp.resize(m_iIMMState, 2);
	gTmp << m_iPeriod / 2, 0,
		1, 0,
		0, m_iPeriod / 2,
		0, 1,
		0, 0,
		0, 0;
	immg.push_back(gTmp);
	gTmp << m_iPeriod*m_iPeriod / 4, 0,
		m_iPeriod / 2, 0,
		0, m_iPeriod*m_iPeriod / 4,
		0, m_iPeriod / 2,
		1, 0,
		0, 1;
	immg.push_back(gTmp);
	immg.push_back(gTmp);

	MatrixXd qTmp;
	qTmp.resize(2, 2);
	qTmp.setZero();
	immq.push_back(qTmp);
	qTmp.setIdentity()*0.001;
	immq.push_back(qTmp);
	qTmp.setIdentity()*0.014;
	immq.push_back(qTmp);

	immH.resize(2, m_iIMMState);
	immH << 1, 0, 0, 0, 0, 0,
		0, 0, 1, 0, 0, 0;

	immR.resize(2, 2);
	immR.setIdentity();
	immR *= RDelta*RDelta;

	immCmean.resize(1, m_iModelNum);
	immCmean.setZero();

	for (int i = 0; i < m_iModelNum; i++){
		immCmean = immCmean + immP.row(i)*immU(i);
	}

	immMu.resize(m_iModelNum, m_iModelNum);
	immMu.setZero();
	//immU(i) ��һά�ȣ��Զ���ά
	for (int i = 0; i < m_iModelNum; i++){
		immMu.row(i) = immU(i) * immP.row(i).array() / immCmean.array();
	}
}


void imm::test()
{
	int montoCarloNum = 5;
	int totalPoint = 900 / m_iPeriod;

	Eigen::MatrixXd zxy, zx, zy, zxEst, zyEst;
	zxy.resize(2 * totalPoint, montoCarloNum);
	zx.resize(totalPoint*montoCarloNum, 1);
	zx.setZero();
	zy.resize(totalPoint*montoCarloNum, 1);
	zy.setZero();

	zxEst.resize(totalPoint, montoCarloNum);
	zxEst.setZero();
	zyEst.resize(totalPoint, montoCarloNum);
	zyEst.setZero();

	//����matlab�����������ļ�
	readTestData("./dataset/immTestData.txt", zxy);
	//cout << zxy.topRows(50) << endl;
	zxy.resize(totalPoint, 2 * montoCarloNum);
	//cout << zxy.topRows(10) << endl;
	zx = zxy.block(0, 0, totalPoint, montoCarloNum);
	zy = zxy.block(0, montoCarloNum, totalPoint, montoCarloNum);
	//cout << zx.topRows(5) << endl;
	//cout << zy.topRows(5) << endl;

	for (int i = 0; i < montoCarloNum; i++){
		zxEst(0, i) = zx(0, i);
		zyEst(0, i) = zy(0, i);
		zxEst(1, i) = zx(1, i);
		zyEst(1, i) = zy(1, i);
		double vx = (zx(1, i) - zx(0, i)) / m_iPeriod;
		double vy = (zy(1, i) - zy(0, i)) / m_iPeriod;

		for (int j = 2; j < totalPoint; j++){//�����в�����ѭ��
			//1�����ù���ֵz_est.(��һ֡��ֵ)
			setEstimateData(zx(j-1, i), vx, zy(j-1, i), vy);
			//2��������ֵz
			setMeasureData(zx(j, i), zy(j, i));
			if (j == 2){
				//3�������ÿ������˲��õ�����ģ�͵Ĺ���ֵ��Э�������
				updateKF();
				zxEst(j, i) = x_est(0, 0);
				zyEst(j, i) = x_est(2, 0);
				//cout << zxEst(j, i) << endl;
			}
			else {
				if (j == 3){
					imm_x_est_exd.block(0, 0, 4, 1) = x_est;
					imm_p_est_exd.block(0, 0, 4, 4) = p_est;
					//4��������ģ�͸���ֵ
					for (int k = 0; k < m_iModelNum; k++){
						imm_xn_est.at(k) = imm_x_est_exd;
						imm_pn_est.at(k) = imm_p_est_exd;
					}
				}
				//5��imm�˲�����
				//[x_est, p_est, xn_est, pn_est, u] = imm(xn_est, pn_est, T, z, Delta, u);
				updateIMM();
				zxEst(j, i) = imm_x_est_exd(0, 0);
				zyEst(j, i) = imm_x_est_exd(2, 0);
			}
		}
	}

	//������õ�������д���ļ��У���matlab����
	writeResultData(zxEst.col(montoCarloNum - 1), "./dataset/immEstDatax.txt");
	writeResultData(zyEst.col(montoCarloNum - 1), "./dataset/immEstDatay.txt");

	//������
	MatrixXd ex, ey;
	ex.resize(totalPoint, 1);
	ey.resize(totalPoint, 1);
	for (int i = 0; i < totalPoint; i++){
		ex(i, 0) = (zx.row(i) - zxEst.row(i)).sum() / montoCarloNum;
		ey(i, 0) = (zy.row(i) - zyEst.row(i)).sum() / montoCarloNum;
	}
}

void imm::updateKF(){
	cout << x_est << endl;
	//cout << phi << endl;
	MatrixXd x_pre = phi*x_est;
	cout << x_pre << endl;
	MatrixXd p_pre = phi*p_est*phi.transpose();
	cout << p_pre << endl;
	MatrixXd kg = p_pre*h.transpose()*(h*p_pre*h.transpose() + R).inverse();
	cout << kg << endl;
	x_est = x_pre + kg*(z - h*x_pre);
	cout << x_est << endl;
	cout << endl;
	MatrixXd I = MatrixXd::Identity(4, 4);
	p_est = (I - kg*h)*p_pre;
	//cout << p_est << endl;

}

void imm::updateIMM()
{
	//1���������룬���ݵ�ǰimm������xn��pn����imm�ڲ��������ת�ƾ��󣬼��㵱ǰÿ��ģ��xn��pn
	vector<MatrixXd> imm_x0, imm_p0;
	//����ÿ��ģ���µĳ�ʼֵ�ͳ�ʼ���ʾ���
	for (int j = 0; j < m_iModelNum; j++){
		MatrixXd tmpx0;
		tmpx0.resize(m_iIMMState, 1);
		tmpx0.setZero();
		MatrixXd tmpp0;
		tmpp0.resize(m_iIMMState, m_iIMMState);
		tmpp0.setZero();
		
		//���㵱ǰģ����ÿ��immMu����ת�ƾ����µĹ���ֵxn��
		for (int i = 0; i < m_iModelNum; i++){
			tmpx0 = tmpx0 + imm_xn_est.at(i)*immMu(i,j);
		}
		imm_x0.push_back(tmpx0);//x0�൱�ڵ�ǰģ���µ�ÿ����ֵ��
		//cout << tmpx0 << endl;
		for (int i = 0; i < m_iModelNum; i++){
			tmpp0 = tmpp0 + immMu(i, j)*(imm_pn_est.at(i)
				+ (imm_xn_est.at(i) - imm_x0.at(j)) 
				* (imm_xn_est.at(i) - imm_x0.at(j)).transpose());
		}
		imm_p0.push_back(tmpp0);
		//cout << tmpp0 << endl;
	}

	//2��ģ�������˲�
	vector<MatrixXd> imm_x_pre, imm_p_pre, imm_kg;
	for (int j = 0; j < m_iModelNum; j++){
		//cout << immPhi.at(j) << endl;
		//cout << imm_x0.at(j) << endl;
		MatrixXd x_pre_tmp = immPhi.at(j)*imm_x0.at(j);
		MatrixXd immQ = immg.at(j)*immq.at(j)*immg.at(j).transpose();
		MatrixXd p_pre_tmp = immPhi.at(j)*imm_p0.at(j)*immPhi.at(j).transpose() + immQ;
		//cout << p_pre_tmp << endl;
		//cout << immR << endl;
		//cout << immH << endl;
		//cout << (immH*p_pre_tmp*immH.transpose() + immR).inverse() << endl;
		MatrixXd kg_tmp = p_pre_tmp*immH.transpose()
			*(immH*p_pre_tmp*immH.transpose() + immR).inverse();
		imm_xn_est.at(j) = x_pre_tmp + kg_tmp * (z - immH * x_pre_tmp);
		//cout << kg_tmp << endl;
		//cout << x_pre_tmp << endl;
		//cout << imm_xn_est.at(j) << endl;
		//cout << p_pre_tmp << endl;
		
		MatrixXd I = MatrixXd(m_iIMMState, m_iIMMState).setIdentity();
		imm_pn_est.at(j) = (I - kg_tmp * immH) * p_pre_tmp;
		//cout << imm_pn_est.at(j) << endl;
		imm_x_pre.push_back(x_pre_tmp);
		imm_p_pre.push_back(p_pre_tmp);
		imm_kg.push_back(kg_tmp);
	}

	//3��ģ�͸��ʸ���
	MatrixXd immA;//
	immA.resize(1, m_iModelNum);//��Ȼ��������ģ��ֵ
	immA.setZero();
	for (int j = 0; j < m_iModelNum; j++){
		//cout << z << endl;
		//cout << imm_x_pre.at(j) << endl;
		MatrixXd immV_tmp = z - immH * imm_x_pre.at(j);//��Ϣֵ
		MatrixXd immS_tmp = immH*imm_p_pre.at(j)*immH.transpose() + immR;//�۲�Э�������
		double immN = immS_tmp.cols()/ 2;
		double coeff = 1 / (pow(2 * gc_dPi, immN) * sqrt(immS_tmp.determinant()));
		//cout << coeff << endl;
		MatrixXd value = coeff * (-0.5*immV_tmp.transpose()*immS_tmp.inverse()*immV_tmp).array().exp();
		//cout << immV_tmp << endl;
		//cout << immS_tmp << endl;
		//cout << value << endl;
		immA(0, j) = value(0, 0);
	}
	//��һ��
	MatrixXd immc = immA.array() * immCmean.array().sum();
	immU = immA.array()*immCmean.array() / immc.array();
	//cout << immU << endl;

	//4���������
	MatrixXd immxn, immpn;
	immxn.resize(m_iIMMState, 1);
	immxn.setZero();
	immpn.resize(m_iIMMState, m_iIMMState);
	immpn.setZero();
	for (int j = 0; j < m_iModelNum; j++){
		immxn = immxn + imm_xn_est.at(j) * immU(j);
		//cout << immxn << endl;
	}
	for (int j = 0; j < m_iModelNum; j++){
		immpn = immpn + immU(j) * (imm_pn_est.at(j)
			+ (imm_xn_est.at(j) - immxn)*(imm_xn_est.at(j) - immxn).transpose());
	}
	imm_x_est_exd = immxn;
	imm_p_est_exd = immpn;
}


void imm::readTestData(string filename, Eigen::MatrixXd& data)
{
	data.resize(data.rows() *data.cols() , 1);
	fstream fin;
	vector<double> testData;
	fin.open(filename, ios::in);
	if (!fin){
		cout << "open file failed" << endl;
		return;
	}
	else{
		int i = 0;
		while (!fin.eof()){
			fin >> data(i,0);
			i++;
			if (i >= data.rows() * data.cols()){
				break;
			}
		}
	}
	data.resize(data.rows(), data.cols());
	fin.close();
}

//����д������
void imm::writeResultData(Eigen::MatrixXd data, string filename)
{
	ofstream fout;
	fout.open(filename, ios::out);
	fout.setf(ios::fixed, ios::floatfield);
	fout.precision(4);

	for (int i = 0; i < data.cols(); i++){
		for (int j = 0; j < data.rows(); j++){
			fout << data(j,i) << endl;
		}
	}

	fout.close();
}
}