#include "JPDA.h"

namespace JPDA{
void JPDA::init()
{
	//��ʼ
	initFlag = 0;
	//����
	m_dPeriod = 1.0;
	//������
	pd = 1.0;
	//��ȷ������������ŵĸ���
	pg = 0.99;
	//����
	gSigma = 9.21;
	//ÿ����λ�����km^2���ڲ���lambda���Ӳ�
	lambda = 2;
	gamma = lambda * exp(-6);

	targetDelta.resize(1, 2);
	targetDelta << 100, 100;

	phi.resize(4, 4);
	phi << 1, m_dPeriod, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, m_dPeriod,
		0, 0, 0, 1;
	h.resize(2, 4);
	h << 1, 0, 0, 0,
		0, 0, 1, 0;
	R.resize(2, 2);
	R << targetDelta(1)*targetDelta(1), 0,
		0, targetDelta(1)*targetDelta(1);
	Q.resize(2, 2);
	Q << 4, 0,
		0, 4;
	G.resize(4, 2);
	G << m_dPeriod*m_dPeriod / 2, 0,
		m_dPeriod, 0,
		0, m_dPeriod*m_dPeriod / 2,
		0, m_dPeriod;

	//����ʵ�Ļ����У��Ǻ����ĸ���
	m_iTargetNum = 2;

	//��ʼ��Э�������
	MatrixXd tmpP;
	tmpP.resize(4, 4);
	tmpP << targetDelta(1)*targetDelta(1), 0, 0, 0,
		0, 0.01, 0, 0,
		0, 0, targetDelta(1)*targetDelta(1), 0,
		0, 0, 0, 0.01;
	for (int i = 0; i < m_iTargetNum; i++){
		p.push_back(tmpP);
	}
	//��ʼ������ֵ��ʼֵ
	for (int i = 0; i < m_iTargetNum; i++){
		MatrixXd tmp;
		tmp.resize(4, 1);
		tmp.setZero();
		xEstimate.push_back(tmp);

		MatrixXd xPredictTmp;
		xPredictTmp.resize(4, 1);
		xPredictTmp.setZero();
		xPredict.push_back(xPredictTmp);

		MatrixXd pPredictTmp;
		pPredictTmp.resize(4, 4);
		pPredictTmp.setZero();
		pPredict.push_back(pPredictTmp);

		MatrixXd zPredictTmp;
		zPredictTmp.resize(2, 1);
		zPredictTmp.setZero();
		zPredict.push_back(zPredictTmp);

		MatrixXd STmp;
		STmp.resize(2, 2);
		STmp.setZero();
		S.push_back(STmp);

		ellipseVol.push_back(0);

	}


}

void JPDA::test()
{
	//tmpData��Ϊÿ����ÿ���������۲�����ŷ�Χ�ڵ�Ŀ���������������Ӳ���
	vector<MatrixXd> tmpData;
	readTxtStreamData("./dataset/jpdaTestData.txt", tmpData);
	MatrixXd initObs;
	initObs.resize(4, 2);
	initObs << 1500, 500,
			   300, 400, 
			   500, 1500,
			   400, 300;

	//����ÿ��ʱ�̵�
	MatrixXd writeData;
	writeData.resize(4, m_iTargetNum*tmpData.size());
	writeData.setZero();
	for (int t = 0; t < tmpData.size(); t++){
		//��ʼ������
		if (t == 0){
			setEstimate(initObs);
			initFlag = 1;
		}
		//ִ��jpda�Ĳ���
		update(tmpData.at(t));
		initFlag = 0;
		//����Ϊִ��д�����ݵ��ļ�����
		vector<MatrixXd> data;
		getEstimate(data);
		for (int i = 0; i < data.size(); i++){
			writeData.col(m_iTargetNum*t + i) = data.at(i);
			//cout << writeData.col(i) << endl;
		}
	}
	writeTxtStreamData(writeData, "./dataset/jpdaEstData.txt");
}

void JPDA::update(MatrixXd measureData)
{
	
	//1����ʼ��Э�������
	for (int i = 0; i < m_iTargetNum; i++){
		MatrixXd xPredictTmp;
		if (initFlag == 0){
			xPredictTmp = phi*xEstimate.at(i);
		}
		else{
			xPredictTmp = xEstimate.at(i);
		}
		//cout << xPredictTmp << endl;
		xPredict.at(i) = xPredictTmp;
		MatrixXd pPredictTmp = phi*p.at(i)*phi.transpose() + G*Q*G.transpose();
		//cout << pPredictTmp << endl;
		pPredict.at(i) = pPredictTmp;
		MatrixXd zPredictTmp = h*xPredict.at(i);
		//cout << zPredictTmp << endl;
		zPredict.at(i) = zPredictTmp;
		R << targetDelta(i)*targetDelta(i), 0,
			0, targetDelta(i)*targetDelta(i);
		MatrixXd SPredictTmp = h*pPredict.at(i)*h.transpose() + R;
		//cout << SPredictTmp << endl;
		S.at(i) = SPredictTmp;
		double ellipseVolTmp = gc_dPi*gSigma*sqrt(S.at(i).determinant());
		//cout << ellipseVolTmp << endl;
		//cout << endl;
		ellipseVol.at(i) = ellipseVolTmp;
	}


	//1������۲�ȷ�Ͼ���Q2����Ϊ���㵱ǰ���ŷ�Χ�ڵĹ۲�Ŀ��ʹ��ں����Ĺ�������
	int measureNum = measureData.cols();
	vector<VectorXi> q2;
	vector<MatrixXd> y;
	for (int j = 0; j < measureNum; j++){
		int flag = 0;//������ʶ��ǰ�ز��Ƿ��ڲ��ŷ�Χ�ڣ����ڱ��浽y������
		VectorXi tmp;
		tmp.resize(m_iTargetNum + 1);//���һ�б����Ӳ���š�
		tmp.setZero();
		tmp.col(0).setOnes();
		for (int i = 0; i < m_iTargetNum; i++){
			MatrixXd d = measureData.col(j) - zPredict.at(i);
			MatrixXd D = d.transpose() * S.at(i).inverse() * d;
			if (D.determinant() <= gSigma){
				flag = 1; 
				tmp(0) = 1;//��һ��ȫ��Ϊ�Ӳ�
				tmp(i+1) = 1;
			}
		}
		//�����벨�ŷ�Χ��ÿ���ز�����������ȥ�������õĻز�
		if (flag == 1){
			q2.push_back(tmp);
			//cout << tmp << endl;
			y.push_back(measureData.col(j));
			//cout << measureData.col(j) << endl;
		}
	}
	//2�����㻥������aMatrix��������ٷ�������������
	int obsMatNum = q2.size();
	vector<MatrixXd> aMatrix;
	MatrixXd aTmp;
	aTmp.resize(obsMatNum, m_iTargetNum + 1);
	aTmp.setZero();
	aTmp.col(0).setOnes();

	//���ÿ��target�Ļ��������������ٷ�����ֻ����������������¡�
	for (int i = 0; i < obsMatNum; i++){
		aTmp.setZero();
		aTmp.col(0).setOnes();
		if (q2.at(i)(1) == 1){//k+1��ʼ����Ϊ����һ��ȫ���Ӳ��У�����Ҫ�жϡ�
			aTmp(i, 1) = 1;
			aTmp(i, 0) = 0;//�������ˣ��ͽ��Ӳ��У�ƥ�����Ϊ0��
			aMatrix.push_back(aTmp);

			for (int j = 0; j < obsMatNum; j++){
				//���³�ʼ��aTmp.
				aTmp.setZero();
				aTmp.col(0).setOnes();
				if ((i != j) && (q2.at(j)(2) == 1)){
					aTmp(j, 2) = 1;
					aTmp(j, 0) = 0;
					aMatrix.push_back(aTmp);
				}
			}
		}
	}

	for (int i = 0; i < obsMatNum; i++){
		if (q2.at(i)(2) == 1){
			//���³�ʼ��aTmp.
			aTmp.setZero();
			aTmp.col(0).setOnes();
			aTmp(i, 2) = 1;
			aTmp(i, 0) = 0;
			aMatrix.push_back(aTmp);
		}
	}
	//���³�ʼ��aTmp.
	aTmp.setZero();
	aTmp.col(0).setOnes();	
	aMatrix.push_back(aTmp);

	for (int i = 0; i < aMatrix.size(); i++){
		//cout << aMatrix.at(i) << endl;
		//cout << endl;
	}

	//3������������Pr
	int aMatrixNum = aMatrix.size();
	MatrixXd pr;
	pr.resize(1, aMatrixNum);
	pr.setZero();
	for (int i = 0; i < aMatrixNum; i++){
		int falseNum = obsMatNum;
		double n = 1.0;
		//cout << aMatrix.at(i) << endl;
		for (int j = 0; j < obsMatNum; j++){
			int meaIndicator = aMatrix.at(i).row(j).block(0, 1, 1, m_iTargetNum).sum(); //�ο�������ʽ4 - 48
			if (meaIndicator == 1){
				falseNum -= 1;
				for (int k = 0; k < m_iTargetNum ; k++){
					if (aMatrix.at(i)(j, k+1) == 1){
						MatrixXd tmp = y.at(j) - zPredict.at(k);
						MatrixXd b = tmp.transpose() * S.at(k).inverse() * tmp;//b ��һ�� 1*1 ����
						double tmp1 = sqrt((2 * gc_dPi*S.at(k)).determinant());
						double tmp2 = exp(-0.5 * b.determinant());
						n = n / tmp1 * tmp2;
						cout << S.at(k) << endl;
					}
				}
			}
			//cout << n << endl;
		}
		//���������
		double a = 1.0;
		if (pd != 1){
			for (int j = 0; j < m_iTargetNum; j++){
				double targetIndicator = aMatrix.at(i).col(j + 1).sum();
				a = a * pow(pd, targetIndicator) * pow(1 - pd, 1 - targetIndicator);
			}
		}
		double sumEllipseVol = accumulate(ellipseVol.begin(), ellipseVol.end(), 0);
		double a1 = 1.0;
		for (int j = 0; j < falseNum; j++){
			a1 *= j+1;
		}
		pr(0,i) = n*a*a1 / (pow(sumEllipseVol, falseNum));
		//cout << pr(0, i) << endl;
		int aa = 1;
	}
	pr = pr.array() / pr.sum();
	//cout << pr << endl;

	//4�������������U
	MatrixXd u;
	u.resize(obsMatNum + 1, m_iTargetNum);//�����һ��Ϊ�Ӳ��������ݡ�
	u.setZero();
	for (int i = 0; i < m_iTargetNum; i++){
		for (int j = 0; j < obsMatNum; j++){
			for (int k = 0; k < aMatrixNum; k++){
				u(j, i) += pr(k) * aMatrix.at(k)(j, i + 1);
			}
		}
	}
	//��������Ŀ��T�����Ĺ������ʴ���U��m1+1,:),��һ��
	for (int i = 0; i < m_iTargetNum; i++){
		u(obsMatNum, i) = 1 - u.col(i).sum();
	}
	//cout << u << endl;
	//5��kalman�˲�
	vector<MatrixXd> kg;
	for (int i = 0; i < m_iTargetNum; i++){
		pPredict.at(i) = phi * p.at(i) * phi.transpose() + G * Q * G.transpose();
		//cout << pPredict.at(i) << endl << endl;
		MatrixXd kgTmp = pPredict.at(i) * h.transpose() * S.at(i).inverse();
		//cout << kgTmp << endl << endl;
		kg.push_back(kgTmp);
		p.at(i) = pPredict.at(i) - (1 - u(obsMatNum, i))* kgTmp * S.at(i)*kgTmp.transpose();
		//cout << p.at(i) << endl << endl;
	}
	for (int i = 0; i < m_iTargetNum; i++){
		//5.1���¹���ֵ
		MatrixXd xPredictTmp;
		xPredictTmp.resize(4, 1);
		xPredictTmp.setZero();
		for (int j = 0; j < obsMatNum; j++){
			xPredictTmp += u(j, i)*(xPredict.at(i) + kg.at(i) * (y.at(j) - zPredict.at(i)));
		}
		xPredictTmp += u(obsMatNum, i) * xPredict.at(i);
		xEstimate.at(i) = xPredictTmp;//��Ϊ����Ĺ���ֵ
		//cout << xEstimate.at(i) << endl << endl;

		//5.2����Э�������
		MatrixXd tmpB;
		tmpB.resize(4, 4);
		tmpB.setZero();
		for (int j = 0; j < obsMatNum + 1; j++){
			MatrixXd tmp;
			if (j == obsMatNum){
				tmp = xPredict.at(i);
			}
			else{
				tmp = xPredict.at(i) + kg.at(i) * (y.at(j) - zPredict.at(i));
			}
			tmpB += u(j, i)*(tmp*tmp.transpose() - xPredictTmp * xPredictTmp.transpose());
		}
		p.at(i) += tmpB;
		//cout << p.at(i) << endl << endl;
	}


}


//�ַ����ָ��
std::vector<std::string> split(std::string str, std::string pattern)
{
	std::string::size_type pos;
	std::vector<std::string> result;
	str += pattern;//��չ�ַ����Է������
	int size = str.size();

	for (int i = 0; i<size; i++)
	{
		pos = str.find(pattern, i);
		if (pos<size)
		{
			std::string s = str.substr(i, pos - i);
			result.push_back(s);
			i = pos + pattern.size() - 1;
		}
	}
	return result;
}

//���ݸ�ʽ���Զ�����Ϊ���ݷָ�split��һ��Ϊһ��Ŀ���λ��
//��һ������Ϊ��һ��Ŀ�����������λ��Ϊ���з�ʽx,y
void JPDA::readTxtStreamData(string filename, vector<MatrixXd>& data)
{
	ifstream fin;
	fin.open(filename, ios::in);
	if (!fin){
		cout << "open file failed" << endl;
		return;
	}
	else{
		string str;
		vector<string> vecStr;
		while (getline(fin, str, '\n')){
			vecStr.push_back(str);
			//cout << str << endl;
		}
		 
		for (int i = 0; i < vecStr.size(); i++){
			MatrixXd tmp;
			int num;
			vector<string> dataStr = split(vecStr.at(i), ",");
			for (int k = 0; k < dataStr.size()-1; k++){
				
				if (k == 0){
					num = atoi(dataStr.at(0).c_str());
				}
				else{
					tmp.resize(1, 2 * num);
					double tmpData = atof(dataStr.at(k).c_str());
					tmp(0, k-1) = tmpData;
				}
			}
			tmp.resize(2, num);
			data.push_back(tmp);
		}
	}


}

//����д������
void JPDA::writeTxtStreamData(Eigen::MatrixXd data, string filename)
{
	ofstream fout;
	fout.open(filename, ios::out);
	fout.setf(ios::fixed, ios::floatfield);
	fout.precision(4);

	for (int i = 0; i < data.cols(); i++){ 
		if (i % 2 == 0){
			fout << "2 , ";
		}
		fout << data(0, i) << " , " ;
		fout << data(2, i) << " , " ;
		if (i % 2 == 1){
			fout << endl;
		}
	}
	fout.close();
}
}