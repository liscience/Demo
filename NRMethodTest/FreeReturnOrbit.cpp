#include "FreeReturnOrbit.h"
#include <Eigen/SVD>
#include <Eigen/Core>
#include<iostream>
#include<fstream>
#include<cmath>
#include<iomanip>
using namespace std;

//���캯��
CFreeROInit::CFreeROInit()
{
	;
}
//��������
CFreeROInit::~CFreeROInit()
{
	;
}



//**********************************************************
///NR���������ɷ��ع���������������̺���
///@Author   Li Zeyue
///@Date	 2022.7.30
///@Input    
///			 ����Ŀ��������������� LLOModOrbElem
///          ��Ʊ��� x
///          �Һ���Ŀ������� Y_tar
///          ����� Out_Ylim
///          ��Ʊ����ο����� int_interval
///          ��Ʊ���������Χ Iter_interval
///          ����������� k
///          �������� Iternum
///@Output
///          �Ƿ������Ĳ���ֵ IsCon
//***********************************************************
void CFreeROInit::NRFreeRTOInit(const CModOrbElem& LLOModOrbElem, CVector<double>& x,
	CVector<double> Y_tar, CVector<double> Out_Ylim, CVector<double> int_interval,
	CMatrix<double> Iter_interval, CVector<double> k, int Iternum, bool&IsCon)
{

	double jd1, jd2;                 //������
	double deltaV;                   //�ٶ�����
	CModOrbElem TLI_LTOModElem;      //����ʱ��ת�ƹ������
	double timeLOI2ROI, backwardTime;
	double phaseError, ascendNode;
	CCoord3 LCO_PRLPosMJ2000, LCO_PRLVelMJ2000;
	CCTime tempTime;
	CCoord3 PosEJ2000, VelEJ2000;
	CPropgateFunc findBackwardtimeFunc;
	std::vector<double> timeArray;
	std::vector<CCoord3> posArray;
	std::vector<CCoord3> velArray;

	CCoord3 pos, vel;
	COrbit prop;

	findBackwardtimeFunc.mLaunchSite.m_Lat = mLaunchSiteLat;
	findBackwardtimeFunc.mLaunchSite.m_Lon = mLaunchSiteLon;
	findBackwardtimeFunc.mLaunchSite.m_Rad = mLaunchHeight + AsCEarthRadius;
	findBackwardtimeFunc.mIsDescend = true;
	//������ʼ�ͽ���ʱ��
	AsTimeToJD(mPRLTime[0], jd1);
	AsTimeToJD(mPRLTime[1], jd2);
	int Optxnum = x.GetSize();
	mLLOModElem = LLOModOrbElem;
	for (double i = jd1; i <= jd2; i = i + mStep*AsCSecToDay)
	{
		AsJDToTime(i, mLOITime);
		mLLOModElemNew = mLLOModElem;
		//�����µĹ������
		AsModOrbElemToCart(mLLOModElem, AsCMoonGrav, pos, vel);
		propSelenocentricElem.m_Epoch = i;
		prop.LunarOrbitPropVarStep(propSelenocentricElem, mStep, 1e-13, 1, pos, vel, timeArray, posArray, velArray);
		AsCartToModOrbElem(pos, vel, AsCMoonGrav, mLLOModElem);
		NRMethod(&CFreeROInit::ConJud, &CFreeROInit::SenCalFun, &CFreeROInit::SenCal,
			&CFreeROInit::SenCalFun, x, Y_tar, Out_Ylim,
			int_interval, Iter_interval, k, Iternum, IsCon);
		if (IsCon != 0)
		{
			//�������
			mFreeReturnOrbitInitCalOutput.PRLTime.push_back(mLOITime);
			double deltat;
			deltat = AsTrueToTPP(mLTOTrueA, mLSIModElemEJ2000.m_PeriRad / (1 - mLSIModElemEJ2000.m_Ecc), mLSIModElemEJ2000.m_Ecc, AsCEarthGrav);

			mTLITime.AddSec(deltat);
			mFreeReturnOrbitInitCalOutput.TLILaunchTime.push_back(mTLITime);

			//��¼���뱼�¹�����
			mFreeReturnOrbitInitCalOutput.LTOInc.push_back(mLSIModElemEJ2000.m_I);
			//��¼α����
			mFreeReturnOrbitInitCalOutput.PRLPseudoLat.push_back(mPRLPseudoLat);
			mFreeReturnOrbitInitCalOutput.PRLPseudoLong.push_back(mPRLPseudoLon);
			mFreeReturnOrbitInitCalOutput.PRLEccentric.push_back(mPRLEccentric);
			mFreeReturnOrbitInitCalOutput.PRLVelInc.push_back(mPRLVelInc);
			//��¼���Ĺ������
			mFreeReturnOrbitInitCalOutput.LCO_PRLModElem.push_back(mLCO_PRLModElem);
			mFreeReturnOrbitInitCalOutput.LSIModElem.push_back(mLSIModElemEJ2000);
			mFreeReturnOrbitInitCalOutput.LSOModElem.push_back(mLSOModElemEJ2000);
			mFreeReturnOrbitInitCalOutput.TimeofFlight.push_back(mTofTemp - deltat);

			mFreeReturnOrbitInitCalOutput.VCPH.push_back(mLSOModElemEJ2000.m_PeriRad - AsCEarthRadius);
			mFreeReturnOrbitInitCalOutput.ROIInc.push_back(mLSOModElemEJ2000.m_I);

			TLI_LTOModElem = mLSIModElemEJ2000;
			TLI_LTOModElem.m_TrueA = mLTOTrueA;
			Cal1Impulse(TLI_LTOModElem, 0, deltaV);
			mFreeReturnOrbitInitCalOutput.TLI_LTOModElem.push_back(TLI_LTOModElem);
			mFreeReturnOrbitInitCalOutput.TLIDeltaV.push_back(deltaV);

			AsModOrbElemToCart(mLCO_PRLModElem, AsCMoonGrav, LCO_PRLPosMJ2000, LCO_PRLVelMJ2000);
			//������ս��ص㴦��λ���ٶ�
			CalHPLOI2VCP(mLOITime, LCO_PRLPosMJ2000, LCO_PRLVelMJ2000, tempTime, PosEJ2000, VelEJ2000, timeArray, posArray, velArray);

			//����ROI�ĵ���γ��
			CalROIECFLatLon(timeArray, posArray, velArray, mLOITime, mAtmospH, timeLOI2ROI, mROIECFLat, mROIECFLon);
			tempTime = mLOITime;
			tempTime.AddSec(timeLOI2ROI);
			mFreeReturnOrbitInitCalOutput.ROITime.push_back(tempTime);
			mFreeReturnOrbitInitCalOutput.ROILat.push_back(mROIECFLat);
			mFreeReturnOrbitInitCalOutput.ROILon.push_back(mROIECFLon);

			//���µ��Ķ�
			if (fabs(AsNegPIToPI(mLSIModElemEJ2000.m_ArgPeri)) <= AsCHalfPI)
			{
				mFreeReturnOrbitInitCalOutput.IsLTODesend.push_back(false);//���µ��Ķ�����
			}
			else
			{
				mFreeReturnOrbitInitCalOutput.IsLTODesend.push_back(true);//���µ��Ķν���
			}
			//�������Ķ�
			if (fabs(AsNegPIToPI(mLCO_PRLModElem.m_ArgPeri)) <= AsCHalfPI)
			{
				mFreeReturnOrbitInitCalOutput.IsLDODesend.push_back(false);//���Ķ�����
			}
			else
			{
				mFreeReturnOrbitInitCalOutput.IsLDODesend.push_back(true);//���Ķν���
			}
			//���ص��Ķ�
			if (fabs(AsNegPIToPI(mLSOModElemEJ2000.m_ArgPeri)) <= AsCHalfPI)
			{
				mFreeReturnOrbitInitCalOutput.IsRTODesend.push_back(false);//���Ķ�����
			}
			else
			{
				mFreeReturnOrbitInitCalOutput.IsRTODesend.push_back(true);//���Ķν���
			}
			goto OUTLOOP;
		}
	OUTLOOP:
		continue;
	}
}

//**********************************************************
///NR���������ɷ��ع���߾��ȼ��������̺���
///@Author   Li Zeyue
///@Date	 2022.7.30
///@Input    
///			 ����Ŀ��������������� LLOModOrbElem
///          ��Ʊ��� x
///          �Һ���Ŀ������� Y_tar
///          ����� Out_Ylim
///          ��Ʊ����ο����� int_interval
///          ��Ʊ���������Χ Iter_interval
///          ����������� k
///          �������� Iternum
///@Output
///          �Ƿ������Ĳ���ֵ IsCon
//***********************************************************
void CFreeROInit::NRFreeRTOHP(const CModOrbElem& LLOModOrbElem, CVector<double>& x,
	CVector<double> Y_tar, CVector<double> Out_Ylim, CVector<double> int_interval,
	CMatrix<double> Iter_interval, CVector<double> k, int Iternum, bool&IsCon)
{

	double jd1, jd2;                 //������
									 //������ʼ�ͽ���ʱ��
	AsTimeToJD(mPRLTime[0], jd1);
	AsTimeToJD(mPRLTime[1], jd2);
	int size;
	size = (jd2 - jd1) / (mStep*AsCSecToDay) + 1;
	//���ø߾�������ṹ���С

	mFreeRTONRCalOutput.TLILaunchTime.resize(size);
	mFreeRTONRCalOutput.TLI_LTOModElem.resize(size);
	mFreeRTONRCalOutput.TLI_LTOModElemE.resize(size);
	mFreeRTONRCalOutput.PRL_LCOModElem.resize(size);
	mFreeRTONRCalOutput.LLOModElem.resize(size);
	mFreeRTONRCalOutput.ReentryGamma.resize(size);
	mFreeRTONRCalOutput.TLIDeltaV.resize(size);
	mFreeRTONRCalOutput.PRLTime.resize(size);
	mFreeRTONRCalOutput.RTOInc.resize(size);
	mFreeRTONRCalOutput.VCPH.resize(size);
	mFreeRTONRCalOutput.epsi.resize(size);
	mFreeRTONRCalOutput.ROITime.resize(size);
	mFreeRTONRCalOutput.ROIECFLat.resize(size);
	mFreeRTONRCalOutput.ROIECFLon.resize(size);
	mFreeRTONRCalOutput.ROIPhaseError.resize(size);
	mFreeRTONRCalOutput.perierror.resize(size);
	mFreeRTONRCalOutput.LOIDeltaV.resize(size);
	mFreeRTONRCalOutput.PRL_LIOModElem.resize(size);
	mFreeRTONRCalOutput.PRL_LIOBeforeModElem.resize(size);
	mFreeRTONRCalOutput.LOI1Time.resize(size);
	mFreeRTONRCalOutput.LOI2Time.resize(size);
	mFreeRTONRCalOutput.LOI2DeltaV.resize(size);
	mFreeRTONRCalOutput.LOI2ModElem.resize(size);
	mFreeRTONRCalOutput.LOI2BeforeModElem.resize(size);
	mFreeRTONRCalOutput.LOI3Time.resize(size);
	mFreeRTONRCalOutput.LOI3DeltaV.resize(size);
	mFreeRTONRCalOutput.LOI3ModElem.resize(size);
	mFreeRTONRCalOutput.LOI3BeforeModElem.resize(size);
	mFreeRTONRCalOutput.LDOModElem.resize(size);
	mFreeRTONRCalOutput.TotalImpulse.resize(size);
	mFreeRTONRCalOutput.TLIPhaseError.resize(size);
	mFreeRTONRCalOutput.backwardTime.resize(size);
	mFreeRTONRCalOutput.ReturnLonError.resize(size);

	CCTime       tempTime, tempTime1;
	double       timeLOI2ROI;
	double       deltV;                 //�ƶ��ٶ�
	double       reentryGamma;          //�����

	CCoord3		 pos, vel;
	CCoord3      PRL_LCOPosMJ2000;
	CCoord3      PRL_LCOVelMJ2000;
	CCoord3      PosEJ2000;
	CCoord3      VelEJ2000;

	double       backwardTime;
	double       phaseError, ascendNode;
	double       RAAN, LLOInc;

	CPropgateFunc findBackwardtimeFunc;
	COrbit prop;
	//��������
	std::vector<double> timeArray;
	std::vector<CCoord> posArray;
	std::vector<CCoord> velArray;

	findBackwardtimeFunc.mLaunchSite.m_Lat = mLaunchSiteLat;
	findBackwardtimeFunc.mLaunchSite.m_Lon = mLaunchSiteLon;
	findBackwardtimeFunc.mLaunchSite.m_Rad = mLaunchHeight + AsCEarthRadius;
	findBackwardtimeFunc.mIsDescend = true;

	int Optxnum = x.GetSize();
	mLLOModElem = LLOModOrbElem;
	for (double i = jd1; i <= jd2; i = i + mStep*AsCSecToDay)
	{
		int num;
		num = (i - jd1) / (mStep*AsCSecToDay);
		AsJDToTime(i, mLOITime);
		mLLOModElemNew = mLLOModElem;
		//�����µĹ������
		AsModOrbElemToCart(mLLOModElem, AsCMoonGrav, pos, vel);
		propSelenocentricElem.m_Epoch = i;
		prop.LunarOrbitPropVarStep(propSelenocentricElem, mStep, 1e-13, 1, pos, vel, timeArray, posArray, velArray);
		AsCartToModOrbElem(pos, vel, AsCMoonGrav, mLLOModElem);
		RAAN = mLLOModElem.m_RAAN;
		LLOInc = mLLOModElem.m_I;
		NRMethod(&CFreeROInit::ConJud, &CFreeROInit::SenHPCalFun, &CFreeROInit::SenCal,
			&CFreeROInit::SenCalFun, x, Y_tar, Out_Ylim,
			int_interval, Iter_interval, k, Iternum, IsCon);
		if (IsCon != 0)
		{
			//��¼�Ż���� 
			Cal1Impulse(mTLIModOrbElem, 0, deltV);
			tempTime = mLOITime;
			double deltat;
			//deltat = AsTrueToTPP(mLTOTrueA, mTLIModOrbElem.m_PeriRad / (1 - mTLIModOrbElem.m_Ecc), mTLIModOrbElem.m_Ecc, AsCEarthGrav);
			AsModOrbElemToCart(mLCO_PRLModElem, AsCMoonGrav, PRL_LCOPosMJ2000, PRL_LCOVelMJ2000);
			// �߾���ģ���ɽ��µ�α״̬���������س���ʱ��״̬
			CalHPLOI2TLI(mLOITime, PRL_LCOPosMJ2000, PRL_LCOVelMJ2000, tempTime, PosEJ2000, VelEJ2000, timeArray, posArray, velArray);
			tempTime = mLOITime;

			//���붯̬�ڴ�			
			CModOrbElem	 tempTLI_LTOModeElem;

			int number = timeArray.size() / 2;
			int counter = 0;
			double *p_time = new double[number];
			double *p_trueanal = new double[number];

			//���TLIʱ��
			for (unsigned int i = timeArray.size() - number; i < timeArray.size(); i++)
			{
				//��ʱ��λ��
				p_time[counter] = timeArray[i];
				AsCartToModOrbElem(posArray[i], velArray[i], AsCEarthGrav, tempTLI_LTOModeElem);
				p_trueanal[counter] = AsNegPIToPI(tempTLI_LTOModeElem.m_TrueA);
				counter++;
			}
			//��ֵ��õ���������ʱ��
			AsInterpPol1D(number, p_trueanal, p_time, 7, false, mLTOTrueA, deltat);
			AsInterpPosVel(timeArray, posArray, velArray, 7, false, deltat, PosEJ2000, VelEJ2000);
			tempTime.AddSec(deltat);
			// TLIʱ�� (UTC) 
			mFreeRTONRCalOutput.TLILaunchTime[num] = tempTime;
			// TLIʱ�̣�LTO�Ĺ����������

			//mFreeRTONRCalOutput.TLILaunchTime[1] = tempTime; 
			AsCartToModOrbElem(PosEJ2000, VelEJ2000, AsCEarthGrav, mTLIModOrbElem);
			mFreeRTONRCalOutput.TLI_LTOModElem[num] = mTLIModOrbElem;
			delete[]p_time;
			delete[]p_trueanal;

			mFreeRTONRCalOutput.TLIDeltaV[num] = deltV;
			mFreeRTONRCalOutput.PRLTime[num] = mLOITime;
			mFreeRTONRCalOutput.PRL_LCOModElem[num] = mLCO_PRLModElem;
			mFreeRTONRCalOutput.RTOInc[num] = mVCPModOrbElem.m_I;
			mFreeRTONRCalOutput.VCPH[num] = mVCPModOrbElem.m_PeriRad - AsCEarthRadius;
			CalReenGamma(mVCPModOrbElem.m_PeriRad - AsCEarthRadius, mAtmospH, mVCPModOrbElem.m_Ecc, reentryGamma);

			mepsi = AsSign(RAAN - mLCO_PRLModElem.m_RAAN)*acos(cos(RAAN - mLCO_PRLModElem.m_RAAN)*sin(mLCO_PRLModElem.m_I)*sin(LLOInc)
				+ cos(mLCO_PRLModElem.m_I)*cos(LLOInc));
			double u = AsZeroTo2PI(atan2(sin(RAAN - mLCO_PRLModElem.m_RAAN)*sin(RAAN) / sin(mepsi),
				(cos(mLCO_PRLModElem.m_I)*cos(mepsi) - cos(LLOInc)) / sin(mLCO_PRLModElem.m_I) / sin(mepsi)));
			if (u < mLCO_PRLModElem.m_ArgPeri)
			{
				u += AsCTwoPI;
			}
			mperierror = sin(u - mLCO_PRLModElem.m_ArgPeri);
			//��ȫ��λ
			//mperierror = (asin(sin(mLCO_PRLModElem.m_RAAN - RAAN)*sin(LLOInc) / sin(mepsi)) - mLCO_PRLModElem.m_ArgPeri);

			mFreeRTONRCalOutput.epsi[num] = mepsi;
			mFreeRTONRCalOutput.perierror[num] = asin(mperierror);
			mFreeRTONRCalOutput.ReentryGamma[num] = reentryGamma;

			AsModOrbElemToCart(mLCO_PRLModElem, AsCMoonGrav, PRL_LCOPosMJ2000, PRL_LCOVelMJ2000);
			//����VCP����λ�á��ٶ�
			CalHPLOI2VCP(mLOITime, PRL_LCOPosMJ2000, PRL_LCOVelMJ2000, tempTime1, PosEJ2000, VelEJ2000, timeArray, posArray, velArray);
			if (timeArray.size() == 0)
			{
				//���������
				cout << "ʱ������Ϊ0" << endl;
			}
			//����ROI�ĵ���γ��
			CalROIECFLatLon(timeArray, posArray, velArray, mLOITime, mAtmospH, timeLOI2ROI, mROIECFLat, mROIECFLon);
			tempTime = mLOITime;
			tempTime.AddSec(timeLOI2ROI);
			//������ս��ص�ͷ��ز���
			mFreeRTONRCalOutput.ROITime[num] = tempTime;
			mFreeRTONRCalOutput.ROIECFLat[num] = mROIECFLat;
			mFreeRTONRCalOutput.ROIECFLon[num] = mROIECFLon;

			CCoord3 ROIPos, ReturnPos;
			CLatLonRad ROIECF(mROIECFLat, mROIECFLon, AsCEarthRadius);
			double    angROI2Return;

			//�����������Ϸ��س�γ�ȶ�Ӧ�ĵ�����
			double returnArcb = asin((tan(mReturnSiteECFLat) / tan(mVCPModOrbElem.m_I)));
			double ROIArcb = asin((tan(mROIECFLat) / tan(mVCPModOrbElem.m_I)));

			if (isnan(returnArcb) || isnan(ROIArcb))
			{
				mFreeRTONRCalOutput.ReturnLonError[num] = 0;
				mFreeRTONRCalOutput.ROIPhaseError[num] = 0;
			}
			else
			{

				mFreeRTONRCalOutput.ReturnLonError[num] = AsNegPIToPI((returnArcb - ROIArcb) + mROIECFLon) - mReturnSiteECFLon;

				CLatLonRad ReturnECF(mReturnSiteECFLat, AsNegPIToPI((returnArcb - ROIArcb) + mROIECFLon), AsCEarthRadius);
				AsLLRToCart(ROIECF, ROIPos);
				AsLLRToCart(ReturnECF, ReturnPos);
				if (mIsRTODesend)
				{
					angROI2Return = -AsSign(ReturnPos[2] - ROIPos[2])*AsAngBetween(ROIPos, ReturnPos);
				}
				else
				{
					angROI2Return = AsSign(ReturnPos[2] - ROIPos[2])*AsAngBetween(ROIPos, ReturnPos);
				}
				mFreeRTONRCalOutput.ROIPhaseError[num] = angROI2Return;
			}

			//���������㾭��
			CCoord TLIpos, TLIvel;
			CCoord TLIposECF, TLIvelECF;
			double jdp1, jdp2;
			CModOrbElem TLI_LTOModElem;
			CMatrix<double> J2000ToECF(3, 3);
			TLI_LTOModElem = mTLIModOrbElem;
			//�ع�ϵ�µĸ���
			AsTimeToJD(mFreeRTONRCalOutput.TLILaunchTime[num], jdp1, jdp2);
			AsModOrbElemToCart(TLI_LTOModElem, AsCDX2GM, TLIpos, TLIvel);

			AsJ2000ToECFMtx(jdp1, jdp2, J2000ToECF);
			TLIposECF = J2000ToECF*TLIpos;
			TLIvelECF = J2000ToECF*TLIvel;

			AsCartToModOrbElem(TLIposECF, TLIvelECF, AsCDX2GM, mFreeRTONRCalOutput.TLI_LTOModElemE[num]);
			mFreeRTONRCalOutput.TLI_LTOModElemE[num].m_ArgPeri = mFreeRTONRCalOutput.TLI_LTOModElemE[num].m_ArgPeri*AsCRadToDeg;
			mFreeRTONRCalOutput.TLI_LTOModElemE[num].m_I = mFreeRTONRCalOutput.TLI_LTOModElemE[num].m_I*AsCRadToDeg;
			mFreeRTONRCalOutput.TLI_LTOModElemE[num].m_RAAN = mFreeRTONRCalOutput.TLI_LTOModElemE[num].m_RAAN*AsCRadToDeg;
			mFreeRTONRCalOutput.TLI_LTOModElemE[num].m_TrueA = mFreeRTONRCalOutput.TLI_LTOModElemE[num].m_TrueA*AsCRadToDeg;

			////����������ྭ
			double siteLon, TLIraan;
			//ȷ��TLIʱ�̷��䳡�������������㾭�ȡ�����
			DeterminPhase(mFreeRTONRCalOutput.TLILaunchTime[num], mTLIModOrbElem, findBackwardtimeFunc.mLaunchSite,
				findBackwardtimeFunc.mIsDescend, findBackwardtimeFunc.mflag, siteLon);
			//��ʼ��ʱ�䡢����
			findBackwardtimeFunc.mTLITime = mFreeRTONRCalOutput.TLILaunchTime[num];
			findBackwardtimeFunc.mLTIOrbElem = mTLIModOrbElem;
			TLIraan = mTLIModOrbElem.m_RAAN;
			//������������1��2���ޡ�TLI��������3��4����
			if (findBackwardtimeFunc.mflag <3)
			{
				TLIraan = AsNegPIToPI(TLIraan);
				siteLon = AsNegPIToPI(siteLon);
			}
			if (siteLon - TLIraan < AsCPI && siteLon - TLIraan > 0)
			{
				//������ʱ��,�ⷶΧΪ0~12Сʱ
				AsNLEBisection(findBackwardtimeFunc, -43200, 0, 1e-6, backwardTime);
			}
			else
			{
				//���ǵ������򣬽����޷���
				if (findBackwardtimeFunc.mflag <= 2)
				{
					findBackwardtimeFunc.mflag += 2;
				}
				else
				{
					findBackwardtimeFunc.mflag -= 2;
				}
				//������ʱ��,�ⷶΧΪ12~24Сʱ
				AsNLEBisection(findBackwardtimeFunc, -86400, -43200, 1e-6, backwardTime);
			}

			//��⹲��ʱ�̵���λ��
			GetLaunchParam(mFreeRTONRCalOutput.TLILaunchTime[num], mTLIModOrbElem, findBackwardtimeFunc.mLaunchSite,
				backwardTime, findBackwardtimeFunc.mIsDescend, phaseError, ascendNode);
			double period = AsPeriRadToPeriod(mTLIModOrbElem.m_PeriRad, 0, AsCEarthGrav);
			if (mThrustTime == 0)
			{
				phaseError = AsZeroTo2PI(mTLIModOrbElem.m_ArgPeri + mTLIModOrbElem.m_TrueA)
					- AsZeroTo2PI(phaseError);
			}
			else
			{
				phaseError = AsZeroTo2PI(mTLIModOrbElem.m_ArgPeri + mTLIModOrbElem.m_TrueA)
					- AsZeroTo2PI(phaseError + 2 * AsCPI / period*(-backwardTime - mThrustTime));
			}

			mFreeRTONRCalOutput.TLIPhaseError[num] = AsZeroTo2PI(phaseError);
			mFreeRTONRCalOutput.backwardTime[num] = backwardTime;


		}
	OUTLOOP:
		continue;
	}
}

// ����Eigen�⣬����SVD�ֽ�ķ���������α�棬Ĭ�����erΪ0
Eigen::MatrixXd pinv_eigen_based(Eigen::MatrixXd& origin, const float er = 0) {
	// ����svd�ֽ�
	Eigen::JacobiSVD<Eigen::MatrixXd> svd_holder(origin,
		Eigen::ComputeThinU |
		Eigen::ComputeThinV);
	// ����SVD�ֽ���
	Eigen::MatrixXd U = svd_holder.matrixU();
	Eigen::MatrixXd V = svd_holder.matrixV();
	Eigen::MatrixXd D = svd_holder.singularValues();

	// ����S����
	Eigen::MatrixXd S(V.cols(), U.cols());
	S.setZero();

	for (unsigned int i = 0; i < D.size(); ++i) {

		if (D(i, 0) > er) {
			S(i, i) = 1 / D(i, 0);
		}
		else {
			S(i, i) = 0;
		}
	}

	// pinv_matrix = V * S * U^T
	return V * S * U.transpose();
}


//**********************************************************
///NR���������̺�����֧�ֳ�������͸߾��ȼ���Ļ�ϲ�����
///@Author   Li Zeyue
///@Date	 2022.7.30
///@Input    
///			 �������жϺ���ָ�� Confunc���Ա
///			 �Һ������㺯��ָ�� Calfunc���Ա
///			 ���жȾ�����㺯��ָ�� SenMa���Ա
///          ��Ʊ��� x
///          �Һ���Ŀ������� Y_tar
///          ����� Out_Ylim
///          ��Ʊ����ο����� int_interval
///          ��Ʊ���������Χ Iter_interval
///          ����������� k
///          �������� Iternum
///@Output
///          �Ƿ������Ĳ���ֵ IsCon
//***********************************************************
void CFreeROInit::NRMethod(Confunc f1, Calfunc f2, SenMa f3, Calfunc f4, CVector<double>& x,
	CVector<double> Y_tar, CVector<double> Out_Ylim, CVector<double> int_interval,
	CMatrix<double> Iter_interval, CVector<double> k, int Iternum, bool&IsCon)
{

	CVector<double> Optx;
	int Optxnum = x.GetSize();
	Optx.Resize(Optxnum);
	CVector<double> Optx_Bingo;
	Optx_Bingo.Resize(Optxnum);
	for (int i = 0; i < Optxnum; i++)
		Optx[i] = x[i];//��������ֵ��ֵ��Optx

	int Resnum = Y_tar.GetSize();
	for (int i = 0; i < Optxnum; i++)
	{
		double val_C = (Iter_interval[i][0] + Iter_interval[i][1]) / 2;
		double wid = (Iter_interval[i][1] - Iter_interval[i][0]) / 2;
		if (Optx[i] >= Iter_interval[i][0] && Optx[i] <= Iter_interval[i][1])
			Optx[i] = Optx[i];
		else
			Optx[i] = (Optx[i] - val_C) - wid*(floor((Optx[i] - val_C) / wid)) + val_C;
	}
	CMatrix<double> sen_result(Optxnum, Resnum);
	CVector<double> NRStep(Optxnum);
	bool SuccessJud;
	for (int j = 0; j < Iternum; j++)
	{

		CVector<double> yback(Resnum);
		(this->*f2)(Optx, yback);
		(this->*f1)(yback, Y_tar, Out_Ylim, SuccessJud);
		IsCon = SuccessJud;
		if (SuccessJud != 0)
		{
			for (int i = 0; i < Optxnum; i++)
			{
				Optx_Bingo[i] = Optx[i];
			}
			break;
		}
		else
			(this->*f3)(f2, Optx, yback, int_interval, k, sen_result);
		CMatrix<double> SenInv(Optxnum, Resnum);
		if (Optxnum == Resnum)
		{

			int m = sen_result.GetSizeRow();
			int n = sen_result.GetSizeCol();
			Eigen::MatrixXd A(m,n);
			for (int i = 0; i < m; i++)
			{
				for (int j = 0; j < n; j++)
					A(i, j) = sen_result[i][j];
			}
			Eigen::MatrixXd Sen_Inv(Optxnum, Resnum);
			Sen_Inv = pinv_eigen_based(A);
			for (int i = 0; i < Optxnum; i++)
			{
				for (int j = 0; j < Resnum; j++)
				{
					SenInv[i][j] = Sen_Inv(i,j);
					/*cout << Sen_Inv(i, j) << endl;*/
				}
					
			}
			//SenInv = sen_result.Inv();//�õ���ǰ�����ֵ���������
			//						  //CMatrix<double> mtxAP;
			//						  //CMatrix<double> mtxU;
			//						  //CMatrix<double> mtxV;
			//						  //bool Squa = sen_result.IsSquare();
			//						  //int rank_sen = sen_result.RankGauss();
			//						  //CMatrix<double> SenInv_H;
			//						  //SenInv_H = sen_result.Transpose();


		}
		else
		{
			//sen_result.Resize(3, 4);
			//sen_result[0][0] = 1; sen_result[0][1] = 1; sen_result[0][2] = 2; sen_result[0][3] = 4;
			//sen_result[1][0] = -1; sen_result[1][1] = -1; sen_result[1][2] = -1; sen_result[1][3] = 1;
			//sen_result[2][0] = 3; sen_result[2][1] = 3; sen_result[2][2] = 2; sen_result[2][3] = 4;
			int m = sen_result.GetSizeRow();
			int n = sen_result.GetSizeCol();
			Eigen::MatrixXd A(m,n);
			for (int i = 0; i < m; i++)
			{
				for (int j = 0; j < n; j++)
					A(i, j) = sen_result[i][j];
			}
			Eigen::MatrixXd Sen_Inv(Optxnum, Resnum);
			Sen_Inv = pinv_eigen_based(A);
			for (int i = 0; i < Optxnum; i++)
			{
				for (int j = 0; j < Resnum; j++)
				{
					SenInv[i][j] = Sen_Inv(i, j);
					cout << Sen_Inv(i, j) << endl;
				}
					
			}
			


		}
		NRStep = SenInv.operator*(yback);//�������Ե�������
		Optx = Optx - NRStep;
		for (int i = 0; i < Optxnum; i++)
			//cout << NRStep[i] << endl;

		for (int i = 0; i < Optxnum; i++)
		{
			double val_C = (Iter_interval[i][0] + Iter_interval[i][1]) / 2;
			double wid = (Iter_interval[i][1] - Iter_interval[i][0]) / 2;
			//double val_C = Iter_interval[i][1] - Iter_interval[i][0];
			//double wid = Iter_interval[i][1] - Iter_interval[i][0];
			if (Optx[i] >= Iter_interval[i][0] && Optx[i] <= Iter_interval[i][1])
				Optx[i] = Optx[i];
			//else if(Optx[i]>Iter_interval[i][1])
			//{
			//	//Optx[i] = (Optx[i] - val_C) - wid*(floor((Optx[i] - val_C) / wid)) + val_C;
			//	Optx[i] = Optx[i]  - wid*(floor((Optx[i] - Iter_interval[i][0]) / wid));
			//}
			else
				//Optx[i] = Optx[i] + wid*(floor((Iter_interval[i][1] - Optx[i]) / wid));
				//Optx[i] = Optx[i] + wid*(floor((Iter_interval[i][1] - Optx[i]) / wid)) ;
				Optx[i] = (Optx[i] - val_C) - wid*(floor((Optx[i] - val_C) / wid)) + val_C;
		}
	}

}
//**********************************************************
///NR���������ɷ��ع�������������ж�
///@Author   Li Zeyue
///@Date	 2022.7.30
///@Input    
///			 �Һ���ʵ�ʼ����� Out_Y
///          �Һ���Ŀ������� Y_tar
///          ����� Out_Ylim
///@Output
///          �Ƿ������Ĳ���ֵ IsCon
//***********************************************************
void CFreeROInit::ConJud(CVector<double> Out_Y, CVector<double> Y_tar, CVector<double> Out_Ylim, bool&IsCon)
{
	int sizenum = Out_Y.GetSize();
	CVector<double> Error;
	Error.Resize(sizenum);
	for (int i = 0; i < sizenum; i++)
		Error[i] = abs(Out_Y[i] - Y_tar[i]);
	for (int i = 0; i < sizenum; i++)
	{
		if (Error[i] > Out_Ylim[i])
		{
			IsCon = 0;
			break;
		}
		else
			IsCon = 1;
	}
}

//**********************************************************
///���ڽ��µ�α�������Һ������㺯����˫����������㣩
///@Author   Li Zeyue
///@Date	 2022.7.30
///@Input    
///			 ��Ʊ��� x
///@Output
///          �Һ��������� result
//***********************************************************
void CFreeROInit::SenCalFun(CVector<double>& x, CVector<double>& result)
{
	//����Ŀ��ֵ
	double timeofFlight;           //����ʱ��
	double tLSI2LOI;               //LSI��LOI��ת��ʱ��
	double timeLSO2VCP;            //LSO��VCP��ת��ʱ��
	double timeLSO2ROI;            //LSO��ROI��ת��ʱ��
	double rr, DELTA;              //������ս��ص�Ĳ���
								   //CVector<double> observations;
								   //observations.Resize(12);
	CCoord3 LCO_PRLPosMJ2000, LCO_PRLVelMJ2000;//���µ㴦LCO�����λ���ٶ�
	CCTime LOITime;                //���¼���ʱ��

								   //���Ż��Ĳ������������Ż�Ŀ��Ĺ��ܺ���
	mPRLPseudoLat = x[0];
	mPRLPseudoLon = x[1];
	mPRLEccentric = x[2];
	mPRLVelInc = x[3];
	//mPRLPseudoLat = -0.016996450218444769;
	//mPRLPseudoLon = 3.1173658551551853;
	//mPRLEccentric = 1.5113342651185437;
	//mPRLVelInc =	0.15867919986113038;
	LOITime = mLOITime;
	//�ɽ��µ�α״̬���������س���ʱ��״̬��˫����ģ�ͣ�
	CalInitPseudoPRL2LSI(LOITime, mLCO_PRLHeight, mPRLPseudoLat, mPRLPseudoLon, mPRLEccentric, mPRLVelInc, mLSIModElemEJ2000);
	// ���������Ӱ���򵽽��ص��ʱ��
	CalTimeFromTrueA(mLSIModElemEJ2000, 0, AsCEarthGrav, timeofFlight);
	timeofFlight = fabs(timeofFlight);
	PerilunePseudoState2MJ2000(LOITime, mLCO_PRLHeight, mPRLPseudoLat, mPRLPseudoLon, mPRLEccentric, mPRLVelInc,
		LCO_PRLPosMJ2000, LCO_PRLVelMJ2000);
	//������µ��ƶ�ǰ�Ĺ������
	AsCartToModOrbElem(LCO_PRLPosMJ2000, LCO_PRLVelMJ2000, AsCMoonGrav, mLCO_PRLModElem);
	CalTimeLSI2LOI(mLCO_PRLModElem, tLSI2LOI);
	mTofTemp = timeofFlight + tLSI2LOI;
	//����ʱ��
	mTLITime = mLOITime;
	mTLITime.AddSec(-mTofTemp);
	//�ɽ��µ�α״̬������ⷵ�س���ʱ��״̬��˫����ģ�ͣ�
	CalInitPseudoPRL2LSO(LOITime, mLCO_PRLHeight, mPRLPseudoLat, mPRLPseudoLon,
		mPRLEccentric, mPRLVelInc, mLSOModElemEJ2000);
	//LSO��VCPʱ��
	CalTimeFromTrueA(mLSOModElemEJ2000, AsCTwoPI - 0.001, AsCEarthGrav, timeLSO2VCP);
	timeLSO2VCP = fabs(timeLSO2VCP);

	//LSO��ROIʱ��
	mTrueAnalROI = AsCTwoPI - acos((mLSOModElemEJ2000.m_PeriRad*(1 + mLSOModElemEJ2000.m_Ecc)
		/ (mAtmospH + AsCEarthRadius) - 1) / mLSOModElemEJ2000.m_Ecc);
	if (isnan(mTrueAnalROI))
	{
		mTrueAnalROI = AsCTwoPI - 0.001;
		timeLSO2ROI = timeLSO2VCP;
	}
	else
	{
		CalTimeFromTrueA(mLSOModElemEJ2000, mTrueAnalROI, AsCEarthGrav, timeLSO2ROI);
	}
	mTofTemp = mTofTemp + tLSI2LOI + timeLSO2ROI;

	//�����ս��ص�߶ȣ�һ��Ϊ5~51km
	rr = mAtmospH + AsCEarthRadius;
	DELTA = AsSqr(rr)*(1 - (1 - AsSqr(mLSOModElemEJ2000.m_Ecc))*(AsSqr(tan(mReentryGamma)) + 1));

	if (mLSOModElemEJ2000.m_Ecc < 0.96)
	{
		//δ������
		mVCPH = 50000;
	}
	else
	{
		//�����ս��ص�߶�
		mVCPH = (rr + sqrt(DELTA)) / (1 + mLSOModElemEJ2000.m_Ecc) / (1 + AsSqr(tan(mReentryGamma))) - AsCEarthRadius;
	}
	mTimeError = 1.0 / mTOFLimit*(mTofTemp - mTof);
	mVCPHError = 1.0 / mHErrorLimit*(mLSOModElemEJ2000.m_PeriRad - AsCEarthRadius - mVCPH);
	mLTOPeriHeightError = 1.0 / mHErrorLimit * (mLSIModElemEJ2000.m_PeriRad - mLTOPeriHeight - AsCEarthRadius);
	mLTOIncError = 1.0 / mIncErrorLimit*(mLSIModElemEJ2000.m_I - mLTOInc);
	mRTOIncError = 1.0 / mIncErrorLimit*(mLSOModElemEJ2000.m_I - mRTOInc);
	result[0] = mLTOPeriHeightError;
	result[1] = mVCPHError;
	result[2] = mLTOIncError;
	result[3] = mRTOIncError;
}
//**********************************************************
///���ڽ��µ����������Һ������㺯����˫����������㣩
///@Author   Li Zeyue
///@Date	 2022.7.30
///@Input    
///			 ��Ʊ��� x
///@Output
///          �Һ��������� result
//***********************************************************
void CFreeROInit::SenCalFunElements(CVector<double>& x, CVector<double>& result)
{
	//����Ŀ��ֵ

	double timeofFlight;           //����ʱ��
	double tLSI2LOI;               //LSI��LOI��ת��ʱ��
	double timeLSO2VCP;            //LSO��VCP��ת��ʱ��
	double timeLSO2ROI;            //LSO��ROI��ת��ʱ��
	double rr, DELTA;              //������ս��ص�Ĳ���
								   //CVector<double> observations;
								   //observations.Resize(12);
	CCoord3 LCO_PRLPosMJ2000, LCO_PRLVelMJ2000;//���µ㴦LCO�����λ���ٶ�
	CCTime LOITime;                //���¼���ʱ��

								   //���Ż��Ĳ������������Ż�Ŀ��Ĺ��ܺ���
	PRLPeriluneEccentric = x[0];
	PRLPeriluneInc = x[1];
	PRLPeriluneRAAN = x[2];
	PRLPeriluneLAOP = x[3];
	LOITime = mLOITime;
	mLCO_PRLModElem.m_PeriRad = mLCO_PRLHeight + AsCMoonRadius;
	mLCO_PRLModElem.m_Ecc = PRLPeriluneEccentric;
	mLCO_PRLModElem.m_I = PRLPeriluneInc;
	mLCO_PRLModElem.m_RAAN = PRLPeriluneRAAN;
	mLCO_PRLModElem.m_ArgPeri = PRLPeriluneLAOP;
	mLCO_PRLModElem.m_TrueA = 0;
	//�ɽ��µ����Ķι���������������س���ʱ��״̬��˫����ģ�ͣ�
	CalInitElementsPRL2LSI(LOITime, mLCO_PRLHeight, PRLPeriluneEccentric, PRLPeriluneInc,
		PRLPeriluneRAAN, PRLPeriluneLAOP, mLSIModElemEJ2000);
	// ���������Ӱ���򵽽��ص��ʱ��
	CalTimeFromTrueA(mLSIModElemEJ2000, 0, AsCEarthGrav, timeofFlight);
	timeofFlight = fabs(timeofFlight);
	//������µ��ƶ�ǰ�Ĺ������
	CalTimeLSI2LOI(mLCO_PRLModElem, tLSI2LOI);
	mTofTemp = timeofFlight + tLSI2LOI;
	//����ʱ��
	mTLITime = mLOITime;
	mTLITime.AddSec(-mTofTemp);
	//�����Ķι������������ⷵ�س���ʱ��״̬��˫����ģ�ͣ�
	CalInitElementsPRL2LSO(LOITime, mLCO_PRLHeight, PRLPeriluneEccentric, PRLPeriluneInc,
		PRLPeriluneRAAN, PRLPeriluneLAOP, mLSOModElemEJ2000);
	//LSO��VCPʱ��
	CalTimeFromTrueA(mLSOModElemEJ2000, AsCTwoPI - 0.001, AsCEarthGrav, timeLSO2VCP);
	timeLSO2VCP = fabs(timeLSO2VCP);

	//LSO��ROIʱ��
	mTrueAnalROI = AsCTwoPI - acos((mLSOModElemEJ2000.m_PeriRad*(1 + mLSOModElemEJ2000.m_Ecc)
		/ (mAtmospH + AsCEarthRadius) - 1) / mLSOModElemEJ2000.m_Ecc);
	if (isnan(mTrueAnalROI))
	{
		mTrueAnalROI = AsCTwoPI - 0.001;
		timeLSO2ROI = timeLSO2VCP;
	}
	else
	{
		CalTimeFromTrueA(mLSOModElemEJ2000, mTrueAnalROI, AsCEarthGrav, timeLSO2ROI);
	}
	mTofTemp = mTofTemp + tLSI2LOI + timeLSO2ROI;

	//�����ս��ص�߶ȣ�һ��Ϊ5~51km
	rr = mAtmospH + AsCEarthRadius;
	DELTA = AsSqr(rr)*(1 - (1 - AsSqr(mLSOModElemEJ2000.m_Ecc))*(AsSqr(tan(mReentryGamma)) + 1));

	if (mLSOModElemEJ2000.m_Ecc < 0.96)
	{
		//δ������
		mVCPH = 50000;
	}
	else
	{
		//�����ս��ص�߶�
		mVCPH = (rr + sqrt(DELTA)) / (1 + mLSOModElemEJ2000.m_Ecc) / (1 + AsSqr(tan(mReentryGamma))) - AsCEarthRadius;
	}
	mTimeError = 1.0 / mTOFLimit*(mTofTemp - mTof);
	mVCPHError = 1.0 / mHErrorLimit*(mLSOModElemEJ2000.m_PeriRad - AsCEarthRadius - mVCPH);
	mLTOPeriHeightError = 1.0 / mHErrorLimit * (mLSIModElemEJ2000.m_PeriRad - mLTOPeriHeight - AsCEarthRadius);
	mLTOIncError = 1.0 / mIncErrorLimit*(mLSIModElemEJ2000.m_I - mLTOInc);
	mRTOIncError = 1.0 / mIncErrorLimit*(mLSOModElemEJ2000.m_I - mRTOInc);


	result[0] = mLTOPeriHeightError;
	result[1] = mVCPHError;
	result[2] = mLTOIncError;
	result[3] = mRTOIncError;

}
//**********************************************************
///���ڽ��µ�α�������Һ������㺯�����߾��ȼ��㣩
///@Author   Li Zeyue
///@Date	 2022.7.30
///@Input    
///			 ��Ʊ��� x
///@Output
///          �Һ��������� result
//***********************************************************
void CFreeROInit::SenHPCalFun(CVector<double>& x, CVector<double>& result)
{
	//α����
	double PRLPseudoLatNew;
	double PRLPseudoLongNew;
	double PRLPeriluneVelIncNew;
	double PRLPeriluneEccentricNew;
	double rr, DELTA;              //������ս��ص�Ĳ���
								   //double epsi, perierror;
								   //double reentryGamma;                  //�����
	CCoord3 PRL_LCOPosMJ2000;             //LCO������µ�λ��
	CCoord3 PRL_LCOVelMJ2000;             //LCO������µ��ٶ�
	CCoord3 PosEJ2000, VelEJ2000;
	CCTime LOITime;                       //���µ����ʱ��
	CCTime tempTime;
	//��������
	std::vector<double> timeArray;
	std::vector<CCoord> posArray;
	std::vector<CCoord> velArray;

	/////////////���Ż��Ĳ������������Ż�Ŀ��Ĺ��ܺ���////////////
	PRLPseudoLatNew = x[0];
	PRLPseudoLongNew = x[1];
	PRLPeriluneEccentricNew = x[2];
	PRLPeriluneVelIncNew = x[3];
	LOITime = mLOITime;

	// ��α�������λ���ٶ�
	PerilunePseudoState2MJ2000(LOITime, mLCO_PRLHeight, PRLPseudoLatNew, PRLPseudoLongNew, PRLPeriluneEccentricNew,
		PRLPeriluneVelIncNew, PRL_LCOPosMJ2000, PRL_LCOVelMJ2000);
	AsCartToModOrbElem(PRL_LCOPosMJ2000, PRL_LCOVelMJ2000, AsCMoonGrav, mLCO_PRLModElem);
	// �߾���ģ���ɽ��µ�α״̬���������س���ʱ��״̬
	CalHPLOI2TLI(LOITime, PRL_LCOPosMJ2000, PRL_LCOVelMJ2000, tempTime, PosEJ2000, VelEJ2000, timeArray, posArray, velArray);
	AsCartToModOrbElem(PosEJ2000, VelEJ2000, AsCEarthGrav, mTLIModOrbElem);
	mTof = (tempTime - LOITime).GetTotalSec();
	// �߾���ģ���ɽ��µ�������ⷵ�ؽ���ʱ��״̬
	// ��ȷ���ĳ�����ɷ��ع��ʱ��������ĺ�����������ֹ�����ֵ����Ŀ��ֵ����Ҫ�������µ�ʱ��
	CalHPLOI2VCP(LOITime, PRL_LCOPosMJ2000, PRL_LCOVelMJ2000, tempTime, PosEJ2000, VelEJ2000, timeArray, posArray, velArray);
	AsCartToModOrbElem(PosEJ2000, VelEJ2000, AsCEarthGrav, mVCPModOrbElem);
	//�����ս��ص�߶ȣ�һ��Ϊ5~51km
	rr = mAtmospH + AsCEarthRadius;
	DELTA = AsSqr(rr)*(1 - (1 - AsSqr(mLSOModElemEJ2000.m_Ecc))*(AsSqr(tan(mReentryGamma)) + 1));

	if (mLSOModElemEJ2000.m_Ecc < 0.96)
	{
		//δ������
		mVCPH = 50000;
	}
	else
	{
		//�����ս��ص�߶�
		mVCPH = (rr + sqrt(DELTA)) / (1 + mLSOModElemEJ2000.m_Ecc) / (1 + AsSqr(tan(mReentryGamma))) - AsCEarthRadius;
	}

	//���ݾ���̶���Ȩ����
	mLTOPeriHeightError = 1.0 / mHErrorLimit * (mTLIModOrbElem.m_PeriRad - mLTOPeriHeight - AsCEarthRadius);
	mVCPHError = 1.0 / mHErrorLimit*(mVCPModOrbElem.m_PeriRad - AsCEarthRadius - mVCPH);
	mLTOIncError = 1.0 / mIncErrorLimit*(mTLIModOrbElem.m_I - mLTOInc);
	mRTOIncError = 1.0 / mIncErrorLimit*(mVCPModOrbElem.m_I - mRTOInc);

	result[0] = mLTOPeriHeightError;
	result[1] = mVCPHError;
	result[2] = mLTOIncError;
	result[3] = mRTOIncError;

}
//**********************************************************
///���ڽ��µ����������Һ������㺯�����߾��ȼ��㣩
///@Author   Li Zeyue
///@Date	 2022.9.7
///@Input    
///			 ��Ʊ��� x
///@Output
///          �Һ��������� result
//***********************************************************
void CFreeROInit::SenHPCalFunElement(CVector<double>& x, CVector<double>& result)
{
	//α����
	double PRLPeriluneEccentric;
	double PRLPeriluneInc;
	double PRLPeriluneRAAN;
	double PRLPeriluneLAOP;

	double rr, DELTA;              //������ս��ص�Ĳ���
	//double epsi, perierror;
	//double reentryGamma;                  //�����
	CCoord3 PRL_LCOPosMJ2000;             //LCO������µ�λ��
	CCoord3 PRL_LCOVelMJ2000;             //LCO������µ��ٶ�
	CCoord3 PosEJ2000, VelEJ2000;
	CCTime LOITime;                       //���µ����ʱ��
	CCTime tempTime;
	//��������
	std::vector<double> timeArray;
	std::vector<CCoord> posArray;
	std::vector<CCoord> velArray;

	/////////////���Ż��Ĳ������������Ż�Ŀ��Ĺ��ܺ���////////////
	PRLPeriluneEccentric = x[0];
	PRLPeriluneInc = x[1];
	PRLPeriluneRAAN = x[2];
	PRLPeriluneLAOP = x[3];
	LOITime = mLOITime;
	mLCO_PRLModElem.m_PeriRad = mLCO_PRLHeight + AsCMoonRadius;
	mLCO_PRLModElem.m_Ecc = PRLPeriluneEccentric;
	mLCO_PRLModElem.m_I = PRLPeriluneInc;
	mLCO_PRLModElem.m_RAAN = PRLPeriluneRAAN;
	mLCO_PRLModElem.m_ArgPeri = PRLPeriluneLAOP;
	mLCO_PRLModElem.m_TrueA = 0;
	AsModOrbElemToCart(mLCO_PRLModElem, AsCMoonGrav, PRL_LCOPosMJ2000, PRL_LCOVelMJ2000);
	// �߾���ģ���ɽ��µ�α״̬���������س���ʱ��״̬
	CalHPLOI2TLI(LOITime, PRL_LCOPosMJ2000, PRL_LCOVelMJ2000, tempTime, PosEJ2000, VelEJ2000, timeArray, posArray, velArray);
	AsCartToModOrbElem(PosEJ2000, VelEJ2000, AsCEarthGrav, mTLIModOrbElem);
	mTof = (tempTime - LOITime).GetTotalSec();
	// �߾���ģ���ɽ��µ�������ⷵ�ؽ���ʱ��״̬
	// ��ȷ���ĳ�����ɷ��ع��ʱ��������ĺ�����������ֹ�����ֵ����Ŀ��ֵ����Ҫ�������µ�ʱ��
	CalHPLOI2VCP(LOITime, PRL_LCOPosMJ2000, PRL_LCOVelMJ2000, tempTime, PosEJ2000, VelEJ2000, timeArray, posArray, velArray);
	AsCartToModOrbElem(PosEJ2000, VelEJ2000, AsCEarthGrav, mVCPModOrbElem);
	//�����ս��ص�߶ȣ�һ��Ϊ5~51km
	rr = mAtmospH + AsCEarthRadius;
	DELTA = AsSqr(rr) * (1 - (1 - AsSqr(mLSOModElemEJ2000.m_Ecc)) * (AsSqr(tan(mReentryGamma)) + 1));

	if (mLSOModElemEJ2000.m_Ecc < 0.96)
	{
		//δ������
		mVCPH = 50000;
	}
	else
	{
		//�����ս��ص�߶�
		mVCPH = (rr + sqrt(DELTA)) / (1 + mLSOModElemEJ2000.m_Ecc) / (1 + AsSqr(tan(mReentryGamma))) - AsCEarthRadius;
	}

	//���ݾ���̶���Ȩ����
	mLTOPeriHeightError = 1.0 / mHErrorLimit * (mTLIModOrbElem.m_PeriRad - mLTOPeriHeight - AsCEarthRadius);
	mVCPHError = 1.0 / mHErrorLimit * (mVCPModOrbElem.m_PeriRad - AsCEarthRadius - mVCPH);
	mLTOIncError = 1.0 / mIncErrorLimit * (mTLIModOrbElem.m_I - mLTOInc);
	mRTOIncError = 1.0 / mIncErrorLimit * (mVCPModOrbElem.m_I - mRTOInc);

	result[0] = mLTOPeriHeightError;
	result[1] = mLTOIncError;
	result[2] = mRTOIncError;
	result[3] = mVCPHError;

}
//**********************************************************
///���ɷ��ع������������жȾ���ǰ���ַ�����
///@Author   Li Zeyue
///@Date	 2022.7.30
///@Input    
///			 �Һ������㺯��ָ�� Calfunc���Ա
///          ��Ʊ��� x
///          ��Ʊ����ο����� int_interval
///          ����������� k
///@Output
///          ���������Ծ��� sen_result
///          �Һ��������� result
///@Return	
//***********************************************************
void CFreeROInit::SenCal(Calfunc f, CVector<double>& x, CVector<double>& result,
	CVector<double> int_interval, CVector<double> k, CMatrix<double>& sen_result)
{
	int Xnum = x.GetSize();
	CVector<double> init_x;
	init_x.Resize(Xnum);
	(this->*f)(x, result);
	//for (int k = 0; k<4; k++)
	//	cout << result[k] << "    " << endl;
	int ynum = result.GetSize();

	CVector<double> delta_result;
	delta_result.Resize(ynum);
	CVector<double> init_result;
	init_result.Resize(ynum);
	for (int i = 0; i < ynum; i++)
		init_result[i] = result[i];

	CVector<double> delta(Xnum);
	for (int i = 0; i < Xnum; i++)
		delta[i] = k[i] * int_interval[i];

	for (int i = 0; i < Xnum; i++)
	{
		for (int j = 0; j < ynum; j++)
		{
			init_x = x;
			init_x[i] = x[i] + delta[i];
			(this->*f)(init_x, delta_result);
			//for (int k = 0; k<4; k++)
			//	cout << delta_result[k] << "    " << endl;
			//for (int k = 0; k<4; k++)
			//	cout << delta[k] << "    " << endl;
			sen_result[j][i] = (delta_result[j] - init_result[j]) / delta[i];

		}
	}
}

//***********************************************************
/// ���µ�α״̬ת��������J2000λ���ٶ�
/// @Author	He Boyong   Zhou Wanmeng
/// @Date	2015-01-21  2018.5.12
/// @Input
/// @Param	timePRL						���µ�ʱ�̣�UTCG
/// @Param	PRLHeight					���µ�߶�,m
/// @Param	PRLPerilunePseudoLat		���µ�������˲ʱ����ϵ�з�λ��γ�ȣ�rad
/// @Param	PRLPerilunePseudoLong		���µ�������˲ʱ����ϵ�з�λ�Ǿ��ȣ�rad
/// @Param	PRLPeriluneEccentric  		���µ�LCO���ƫ���ʣ�1.3���ȼ�����µ��ٶȣ�m/s��һ����2.55e3 m/s
/// @Param	PRLPeriluneVelInc			���µ�α��ǣ�������pi��0.20943951��rad = 180��12��deg
/// @Output
/// @Param	posMJ2000				    ����J2000λ�ã�m
/// @Param	velMJ2000				    ����J2000�ٶȣ�m/s
//***********************************************************
void CFreeROInit::PerilunePseudoState2MJ2000(const CCTime& timePRL, const double& PRLHeight, const double& PRLPerilunePseudoLat,
	const double& PRLPerilunePseudoLong, const double& PRLPeriluneEccentric, const double& PRLPeriluneVelInc,
	CCoord3& posMJ2000, CCoord3& velMJ2000)
{
	CCoord3 posMPI, velMPI;				//���Ľ��µ�����ϵ
	CCoord3 posLVLH, velLVLH;			//���Ĺ��˲ʱ����ϵ
	CMatrix<double> LVLHtoMPI(3, 3);    //LVLH��MPI����ϵ
	CMatrix<double> LVLHtoMJ2000(3, 3); //LVLH��MJ2000����ϵ
	double semiMajor, PRLPeriluneVel;   //�볤�ᡢ���µ��ٶ�
	double prlJd;

	//�볤��
	semiMajor = (PRLHeight + AsCMoonRadius) / (1 - PRLPeriluneEccentric);
	//���µ��ٶ�
	PRLPeriluneVel = sqrt(AsCMoonGrav*(2.0 / (PRLHeight + AsCMoonRadius) - 1.0 / semiMajor));
	//���µ�λ�á��ٶ�
	posMPI[0] = PRLHeight + AsCMoonRadius;
	posMPI[1] = 0;
	posMPI[2] = 0;
	velMPI[0] = 0;
	velMPI[1] = PRLPeriluneVel*cos(PRLPeriluneVelInc);
	velMPI[2] = PRLPeriluneVel*sin(PRLPeriluneVelInc);

	//���˲ʱ����ϵMOI�����µ�����ϵMPI��ת������
	AsEulerToMtx(CEuler(PRLPerilunePseudoLong, -PRLPerilunePseudoLat, 0), 321, LVLHtoMPI);
	posLVLH = LVLHtoMPI.Transpose()*posMPI;
	velLVLH = LVLHtoMPI.Transpose()*velMPI;
	//�µ�ת�Ƴ���ʱ�̶�Ӧ������
	AsTimeToJD(timePRL, prlJd);
	//���ɽ��µ�����ϵ�����˲ʱ����ϵ��ת������
	MoonLVLHToJ2000(timePRL, LVLHtoMJ2000);

	//ת��������J2000����ϵ
	posMJ2000 = LVLHtoMJ2000*posLVLH;
	velMJ2000 = LVLHtoMJ2000*velLVLH;
}
//***********************************************************
/// �ɽ��µ�α״̬���������س���ʱ��״̬��˫����ģ�ͣ�
/// @Author	He Boyong   Zhou Wanmeng
/// @Date	2015-01-21  2018.5.12
/// @Input
/// @Param	timePRL						���µ�ʱ�̣�UTCG
/// @Param	PRLHeight					���µ�߶�,m
/// @Param	PRLPerilunePseudoLat		���µ�������˲ʱ����ϵ�з�λ��γ�ȣ�rad
/// @Param	PRLPerilunePseudoLong		���µ�������˲ʱ����ϵ�з�λ�Ǿ��ȣ�rad
/// @Param	PRLPeriluneEccentric  		���µ�LCO���ƫ���ʣ�1.3���ȼ�����µ��ٶȣ�m/s��һ����2.55e3 m/s
/// @Param	PRLPeriluneVelInc			���µ�α��ǣ�������pi��0.20943951��rad = 180��12��deg
/// @Output
/// @Param	LSIModElemEJ2000      ��������Ӱ����ʱ��EJ2000�������������
/// @Return		
//***********************************************************
void CFreeROInit::CalInitPseudoPRL2LSI(const CCTime& timePRL, const double& PRLHeight, const double& PRLPerilunePseudoLat, const double& PRLPerilunePseudoLong,
	const double& PRLPeriluneEccentric, double PRLPeriluneVelInc, CModOrbElem& LSIModElemEJ2000)
{
	CCoord3 posMJ2000, velMJ2000;
	CModOrbElem PRL_LCOModOrbElem;

	//��α״̬ת��Ϊ���µ�MJ2000λ���ٶ�
	PerilunePseudoState2MJ2000(timePRL, PRLHeight, PRLPerilunePseudoLat,
		PRLPerilunePseudoLong, PRLPeriluneEccentric, PRLPeriluneVelInc, posMJ2000, velMJ2000);
	//��λ���ٶȼ�����µ�LCO�����������
	AsCartToModOrbElem(posMJ2000, velMJ2000, AsCMoonGrav, PRL_LCOModOrbElem);
	//��PRL_LCOModOrbElem��������Ӱ������ڵ��EJ2000�������
	CalInitMJ2000PRL2LSI(timePRL, PRL_LCOModOrbElem, PRL_LCOModOrbElem.m_ArgPeri + PRL_LCOModOrbElem.m_TrueA,
		PRL_LCOModOrbElem.m_Ecc, LSIModElemEJ2000);
}
void CFreeROInit::CalInitElementsPRL2LSI(const CCTime& timePRL, const double& PRLHeight, const double& PRLPeriluneEccentric, const double& PRLPeriluneInc,
	const double& PRLPeriluneRAAN, double PRLPeriluneLAOP, CModOrbElem& LSIModElemEJ2000)
{
	CCoord3 posMJ2000, velMJ2000;
	CModOrbElem PRL_LCOModOrbElem;
	//���µ�LCO�����������
	PRL_LCOModOrbElem.m_PeriRad = PRLHeight + AsCMoonRadius;
	PRL_LCOModOrbElem.m_Ecc = PRLPeriluneEccentric;
	PRL_LCOModOrbElem.m_I = PRLPeriluneInc;
	PRL_LCOModOrbElem.m_RAAN = PRLPeriluneRAAN;
	PRL_LCOModOrbElem.m_ArgPeri = PRLPeriluneLAOP;
	PRL_LCOModOrbElem.m_TrueA = 0;
	//��PRL_LCOModOrbElem��������Ӱ������ڵ��EJ2000�������
	CalInitMJ2000PRL2LSI(timePRL, PRL_LCOModOrbElem, PRL_LCOModOrbElem.m_ArgPeri + PRL_LCOModOrbElem.m_TrueA,
		PRL_LCOModOrbElem.m_Ecc, LSIModElemEJ2000);
}
/// �߾���ģ���ɽ��µ�������ⷵ�ؽ���ʱ��״̬
/// @Author	He Boyong   Zhou Wanmeng
/// @Date	2015-04-12  2018.7.3
/// @Input
/// @Param	timePRL		            ���µ�ʱ��,UTCG
/// @Param	PRL_LCOPosMJ2000		����J2000����ϵ˫���߹�����µ�λ��,m
/// @Param	PRL_LCOVelMJ2000		����J2000����ϵ˫���߹�����µ��ٶ�,m/s
/// @Output
/// @Param	time		            ���ص�ʱ��,UTCG
/// @Param	PosEJ2000				����J2000����ϵLTO������ص�λ��,m
/// @Param	VelJ2000				����J2000����ϵLTO������ص��ٶ�,m/s
/// @Param  timeArray               ��ʷ����
/// @Param  posArray
/// @Param  velArray
/// @Return	
//***********************************************************
bool CFreeROInit::CalHPLOI2VCP(const CCTime& timePRL, const CCoord3& PRL_LCOPosMJ2000, const CCoord3& PRL_LCOVelMJ2000,
	CCTime& time, CCoord3& PosEJ2000, CCoord3& VelEJ2000, std::vector<double>& timeArray,
	std::vector<CCoord3>& posArray, std::vector<CCoord3>& velArray)
{
	double PRLJd;                                //PRLʱ��
	CCTime TT;
	CModOrbElem PRL_LCOModOrbElem;               //���µ�����������������
	CCoord3 lunarPosEJ2000, lunarVelEJ2000;      //�����ڵ���J2000����ϵ�Ĳ���
	CCoord3 PRL_LCOPosEJ2000, PRL_LCOVelEJ2000;  //���µ���Ե����EJ2000λ���ٶ�
	bool    isSuccess;

	//�������
	timeArray.resize(0);
	posArray.resize(0);
	velArray.resize(0);

	//����PRLʱ��
	AsUTCToTT(timePRL, TT);
	AsTimeToJD(TT, PRLJd);

	//�����ڵ���J2000����ϵ�Ĳ���
	AsJplPlanetaryEph(PRLJd, AsEJplMoon, AsEJplEarth, lunarPosEJ2000, lunarVelEJ2000);
	//���µ���Ե����EJ2000λ���ٶ�
	PRL_LCOPosEJ2000 = lunarPosEJ2000 + PRL_LCOPosMJ2000;
	PRL_LCOVelEJ2000 = lunarVelEJ2000 + PRL_LCOVelMJ2000;
	//LCO�������������EJ2000
	AsCartToModOrbElem(PRL_LCOPosEJ2000, PRL_LCOVelEJ2000, AsCEarthGrav, PRL_LCOModOrbElem);

	///////////////���û���ֹͣ����////////////////
	CStateTrigger stopperigee;
	stopperigee.m_TripType = 3;
	stopperigee.m_Condition = 0;  ///< 0=Cross Increasing;1=Cross Decreasing;2=Cross Either
	stopperigee.m_RepeatCount = 0;///< requested count, �����ظ�����
	stopperigee.m_TripValue = 0;  ///< ����ֵ(m, rad)
	stopperigee.m_Tolerance = 1e-6;///< �����������(m, rad)
	double duration = 5184000;

	//���ø߾��ȹ�����ƣ��䲽����������ֹ����
	isSuccess = CalHPOPEphemeris(timePRL, false, PRL_LCOModOrbElem, 0, stopperigee, 2, duration, PosEJ2000, VelEJ2000, timeArray, posArray, velArray);
	time = timePRL;
	time.AddSec(duration);

	return isSuccess;
}


//***********************************************************
/// �ɽ��µ�α״̬������ⷵ�س���ʱ��״̬��˫����ģ�ͣ�
/// @Author	He Boyong   Zhou Wanmeng
/// @Date	2015-01-21  2018.5.12
/// @Input
/// @Param	timePRL						���µ�ʱ�̣�UTCG
/// @Param	PRLHeight					���µ�߶�,m
/// @Param	PRLPerilunePseudoLat		���µ�������˲ʱ����ϵ�з�λ��γ�ȣ�rad
/// @Param	PRLPerilunePseudoLong		���µ�������˲ʱ����ϵ�з�λ�Ǿ��ȣ�rad
/// @Param	PRLPeriluneEccentric  		���µ�LCO���ƫ���ʣ�1.3���ȼ�����µ��ٶȣ�m/s��һ����2.55e3 m/s
/// @Param	PRLPeriluneVelInc			���µ�α��ǣ�������pi��0.20943951��rad = 180��12��deg
/// @Output
/// @Param	LSOModElemEJ2000            �ɳ�����Ӱ����ʱ��EJ2000�������������
/// @Return		
//***********************************************************
void CFreeROInit::CalInitPseudoPRL2LSO(const CCTime& timePRL, const double& PRLHeight, const double& PRLPerilunePseudoLat, const double& PRLPerilunePseudoLong,
	const double& PRLPeriluneEccentric, double PRLPeriluneVelInc, CModOrbElem& LSOModElemEJ2000)
{
	CCTime tempTime, TT;							//���µ�ʱ��
	CCoord3 posMJ2000, velMJ2000;				//���µ�λ�á��ٶ�
	CCoord3 posEJ2000, velEJ2000;				//���ڵ�λ�á��ٶ�
	CCoord3 lunarPosEJ2000, lunarVelEJ2000;     //�����λ�á��ٶ�
	double  timeLOI2LSO;						//�ɳ�Ӱ�����ʱ��
	double  jdp;								//LSOʱ�̶�Ӧ������jd,

	CModOrbElem PRL_LCOModOrbElem;

	//��α״̬ת��Ϊ���µ�MJ2000λ���ٶ�
	PerilunePseudoState2MJ2000(timePRL, PRLHeight, PRLPerilunePseudoLat, PRLPerilunePseudoLong, PRLPeriluneEccentric, PRLPeriluneVelInc, posMJ2000, velMJ2000);
	//��λ���ٶȼ�����µ�LCO�����������
	AsCartToModOrbElem(posMJ2000, velMJ2000, AsCMoonGrav, PRL_LCOModOrbElem);

	//�ο�ʱ�̸�ֵ
	tempTime = timePRL;

	//�����µ㵽����Ӱ������ڵ����Ķ���˫���߷���ʱ��
	CalTimeLSI2LOI(PRL_LCOModOrbElem, timeLOI2LSO);

	//���LSO����λ�ú��ٶ�ʸ��
	AsTwoBodyOrbProp(timeLOI2LSO, AsCMoonGrav, posMJ2000, velMJ2000);

	//������J2000����ϵ����Ӱ����߽�LSOת��Ϊ����J2000����ϵ����
	tempTime.AddSec(timeLOI2LSO);

	//LSOʱ�̶�Ӧ������	
	AsUTCToTT(tempTime, TT);
	AsTimeToJD(TT, jdp);

	//JPL���������ڵ�˲ʱ��������ڵ����λ���ٶ�
	AsJplPlanetaryEph(jdp, AsEJplMoon, AsEJplEarth, lunarPosEJ2000, lunarVelEJ2000);

	//����Ӱ������ڵ����J2000����ϵ��״̬
	posEJ2000 = posMJ2000 + lunarPosEJ2000;
	velEJ2000 = velMJ2000 + lunarVelEJ2000;

	//���LSO����EJ2000�����������
	AsCartToModOrbElem(posEJ2000, velEJ2000, AsCEarthGrav, LSOModElemEJ2000);

}
void CFreeROInit::CalInitElementsPRL2LSO(const CCTime& timePRL, const double& PRLHeight, const double& PRLPeriluneEccentric, const double& PRLPeriluneInc,
	const double& PRLPeriluneRAAN, double PRLPeriluneLAOP, CModOrbElem& LSOModElemEJ2000)
{
	CCTime tempTime, TT;							//���µ�ʱ��
	CCoord3 posMJ2000, velMJ2000;				//���µ�λ�á��ٶ�
	CCoord3 posEJ2000, velEJ2000;				//���ڵ�λ�á��ٶ�
	CCoord3 lunarPosEJ2000, lunarVelEJ2000;     //�����λ�á��ٶ�
	double  timeLOI2LSO;						//�ɳ�Ӱ�����ʱ��
	double  jdp;								//LSOʱ�̶�Ӧ������jd,

	CModOrbElem PRL_LCOModOrbElem;
	PRL_LCOModOrbElem.m_PeriRad = PRLHeight + AsCMoonRadius;
	PRL_LCOModOrbElem.m_Ecc = PRLPeriluneEccentric;
	PRL_LCOModOrbElem.m_I = PRLPeriluneInc;
	PRL_LCOModOrbElem.m_RAAN = PRLPeriluneRAAN;
	PRL_LCOModOrbElem.m_ArgPeri = PRLPeriluneLAOP;
	PRL_LCOModOrbElem.m_TrueA = 0;
	//�ο�ʱ�̸�ֵ
	tempTime = timePRL;

	//�����µ㵽����Ӱ������ڵ����Ķ���˫���߷���ʱ��
	CalTimeLSI2LOI(PRL_LCOModOrbElem, timeLOI2LSO);
	//�ɽ��µ�LCO���������������λ���ٶ�
	AsModOrbElemToCart(PRL_LCOModOrbElem, AsCMoonGrav, posMJ2000, velMJ2000);
	//���LSO����λ�ú��ٶ�ʸ��
	AsTwoBodyOrbProp(timeLOI2LSO, AsCMoonGrav, posMJ2000, velMJ2000);

	//������J2000����ϵ����Ӱ����߽�LSOת��Ϊ����J2000����ϵ����
	tempTime.AddSec(timeLOI2LSO);

	//LSOʱ�̶�Ӧ������	
	AsUTCToTT(tempTime, TT);
	AsTimeToJD(TT, jdp);

	//JPL���������ڵ�˲ʱ��������ڵ����λ���ٶ�
	AsJplPlanetaryEph(jdp, AsEJplMoon, AsEJplEarth, lunarPosEJ2000, lunarVelEJ2000);

	//����Ӱ������ڵ����J2000����ϵ��״̬
	posEJ2000 = posMJ2000 + lunarPosEJ2000;
	velEJ2000 = velMJ2000 + lunarVelEJ2000;

	//���LSO����EJ2000�����������
	AsCartToModOrbElem(posEJ2000, velEJ2000, AsCEarthGrav, LSOModElemEJ2000);

}

bool CFreeROInit::HPSaveFreeRTOResult(const std::string& setFilepath)
{
	std::ofstream outfile;
	CCoord3 PosEJ2000, VelEJ2000;
	CCoord3 PosMJ2000, VelMJ2000;
	CCoord3 lunarPosEJ2000, lunarVelEJ2000;
	CCoord3  PRL_LCOPosEJ2000, PRL_LCOVelEJ2000;

	double  jd;
	CCTime  TT;
	CCTime  EndTime;
	std::string launchTimestr, ROItimestr, PRLtimestr;
	std::string timestr1, timestr2, timestr3, timeend;
	outfile.open(setFilepath.c_str());//�����·���й� ��Ҫ��
	if (!outfile)
	{
		cout << "Open Failed!!!" << endl;
		return false;
	}
	outfile.setf(ios::right);
	outfile.precision(8);
	for (unsigned int i = 0; i<mFreeRTONRCalOutput.TLILaunchTime.size(); i++)
	{
		if (!mIsImpulse && mFreeRTONRCalOutput.TLIDeltaV[i] == 0)
		{
			continue;
		}
		if (mIsImpulse && mFreeRTONRCalOutput.TotalImpulse[i] == 0)
		{
			continue;
		}
		AsModOrbElemToCart(mFreeRTONRCalOutput.TLI_LTOModElem[i], AsCEarthGrav, PosEJ2000, VelEJ2000);
		AsModOrbElemToCart(mFreeRTONRCalOutput.PRL_LCOModElem[i], AsCMoonGrav, PosMJ2000, VelMJ2000);
		AsUTCToTT(mFreeRTONRCalOutput.PRLTime[i], TT);
		AsTimeToJD(TT, jd);
		//���µ�����ڵ���J2000����ϵ�Ĳ���
		AsJplPlanetaryEph(jd, AsEJplMoon, AsEJplEarth, lunarPosEJ2000, lunarVelEJ2000);
		//���µ���Ե����EJ2000λ���ٶ�
		PRL_LCOPosEJ2000 = lunarPosEJ2000 + PosMJ2000;
		PRL_LCOVelEJ2000 = lunarVelEJ2000 + VelMJ2000;

		mFreeRTONRCalOutput.TLILaunchTime[i].ToGregStrCh(-6, launchTimestr);
		mFreeRTONRCalOutput.PRLTime[i].ToGregStrCh(-6, PRLtimestr);
		mFreeRTONRCalOutput.ROITime[i].ToGregStrCh(-6, ROItimestr);
		mFreeRTONRCalOutput.LOI1Time[i].ToGregStrCh(-6, timestr1);
		mFreeRTONRCalOutput.LOI2Time[i].ToGregStrCh(-6, timestr2);
		mFreeRTONRCalOutput.LOI3Time[i].ToGregStrCh(-6, timestr3);

		EndTime = mFreeRTONRCalOutput.PRLTime[i];
		EndTime.AddSec(mPRL2RVDTime);
		EndTime.ToGregStrCh(-6, timeend);

		///�߾��Ƚ��
		outfile
			<< "============================" << endl
			<< "���ɷ��ع���ķ���������" << i << endl
			<< "============================" << endl << endl
			<< "����ת��ʱ�䣨sec��:   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.backwardTime[i] << "   " << endl << endl
			<< "����ת�Ƶ��ļнǣ�deg��:   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLIPhaseError[i] * AsCRadToDeg << "   " << endl << endl
			<< "-----------" << endl
			<< "�����״̬" << endl
			<< "-----------" << endl << endl
			<< "��ʼ��Ԫ(UTCG):   " << endl << endl
			<< setw(18) << setprecision(15) << launchTimestr << "   " << endl << endl

			<< "��ʼ��Ԫ(UTCG):   " << endl << endl
			<< setw(18) << setprecision(15) << launchTimestr << "   " << endl << endl
			<< "������� VVLH��m/s��:   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLIDeltaV[i] << "   "
			<< setw(18) << setprecision(15) << 0 << "  "
			<< setw(18) << setprecision(15) << 0 << endl << endl
			<< "EJ2000λ��(m)���ٶ�(m/s):   " << endl << endl
			<< "x:" << "  " << setw(18) << setprecision(15) << PosEJ2000[0] << endl
			<< "y:" << "  " << setw(18) << setprecision(15) << PosEJ2000[1] << endl
			<< "z:" << "  " << setw(18) << setprecision(15) << PosEJ2000[2] << endl
			<< "vx:" << "  " << setw(18) << setprecision(15) << VelEJ2000[0] << endl
			<< "vy:" << "  " << setw(18) << setprecision(15) << VelEJ2000[1] << endl
			<< "vz:" << "  " << setw(18) << setprecision(15) << VelEJ2000[2] << endl << endl
			<< "EJ2000�����������:   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_PeriRad << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_Ecc << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_I*AsCRadToDeg << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_RAAN*AsCRadToDeg << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_ArgPeri*AsCRadToDeg << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_TrueA*AsCRadToDeg << endl << endl
			<< "EJ2000�������:   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_PeriRad / (1 - mFreeRTONRCalOutput.TLI_LTOModElem[i].m_Ecc) << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_Ecc << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_I*AsCRadToDeg << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_RAAN*AsCRadToDeg << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_ArgPeri*AsCRadToDeg << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_TrueA*AsCRadToDeg << endl << endl
			<< "----------" << endl
			<< "���µ�״̬" << endl
			<< "----------" << endl << endl
			<< "��ʼ��Ԫ(UTCG):   " << endl << endl
			<< setw(18) << setprecision(15) << PRLtimestr << endl << endl
			<< "EJ2000λ��(m)���ٶ�(m/s):   " << endl << endl
			<< "x:" << setw(18) << setprecision(15) << PRL_LCOPosEJ2000[0] << endl
			<< "y:" << setw(18) << setprecision(15) << PRL_LCOPosEJ2000[1] << endl
			<< "z:" << setw(18) << setprecision(15) << PRL_LCOPosEJ2000[2] << endl
			<< "vx:" << setw(18) << setprecision(15) << PRL_LCOVelEJ2000[0] << endl
			<< "vy:" << setw(18) << setprecision(15) << PRL_LCOVelEJ2000[1] << endl
			<< "vz:" << setw(18) << setprecision(15) << PRL_LCOVelEJ2000[2] << endl << endl
			<< "MJ2000�����������:   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LCOModElem[i].m_PeriRad << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LCOModElem[i].m_Ecc << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LCOModElem[i].m_I *AsCRadToDeg << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LCOModElem[i].m_RAAN *AsCRadToDeg << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LCOModElem[i].m_ArgPeri *AsCRadToDeg << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LCOModElem[i].m_TrueA *AsCRadToDeg << endl << endl
			<< "----------" << endl
			<< "�����״̬" << endl
			<< "----------" << endl << endl
			<< "���ع�����(deg):   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.RTOInc[i] * AsCRadToDeg << endl << endl
			<< "���ص㾭�Ȳ�(deg):   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.ReturnLonError[i] * AsCRadToDeg << endl << endl
			<< "��ս��ص�߶ȣ�m��:   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.VCPH[i] << endl << endl
			<< "������Ԫ(UTCG):   " << endl << endl
			<< setw(18) << setprecision(15) << ROItimestr << endl << endl
			<< "�����(deg):   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.ReentryGamma[i] * AsCRadToDeg << endl << endl
			<< "��������γ��(deg):   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.ROIECFLat[i] * AsCRadToDeg << endl << endl
			<< "����������(deg):   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.ROIECFLon[i] * AsCRadToDeg << endl << endl
			<< "�������Ľ�(deg):   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.ROIPhaseError[i] * AsCRadToDeg << endl << endl;
		if (!mIsImpulse)
		{
			outfile << endl;
		}
		else
		{
			outfile
				<< "----------" << endl
				<< "��һ�α��" << endl
				<< "----------" << endl << endl
				<< "�����Ԫ��UTCG��:   " << endl << endl
				<< setw(18) << setprecision(15) << timestr1 << endl << endl
				<< "������� MJ2000��m/s��:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOIDeltaV[i][0] << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOIDeltaV[i][1] << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOIDeltaV[i][2] << endl << endl
				<< "ת��ǰMJ2000�����������:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOBeforeModElem[i].m_PeriRad << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOBeforeModElem[i].m_Ecc << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOBeforeModElem[i].m_I* AsCRadToDeg << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOBeforeModElem[i].m_RAAN * AsCRadToDeg << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOBeforeModElem[i].m_ArgPeri* AsCRadToDeg << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOBeforeModElem[i].m_TrueA* AsCRadToDeg << endl << endl
				<< "MJ2000�����������:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOModElem[i].m_PeriRad << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOModElem[i].m_Ecc << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOModElem[i].m_I *AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOModElem[i].m_RAAN *AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOModElem[i].m_ArgPeri *AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOModElem[i].m_TrueA *AsCRadToDeg << endl << endl
				<< "----------" << endl
				<< "�ڶ��α��" << endl
				<< "----------" << endl << endl
				<< "�����Ԫ��UTCG��:   " << endl << endl
				<< setw(18) << setprecision(15) << timestr2 << endl << endl
				<< "������� MJ2000��m/s��:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2DeltaV[i][0] << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2DeltaV[i][1] << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2DeltaV[i][2] << endl << endl
				<< "ת��ǰMJ2000�����������:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2BeforeModElem[i].m_PeriRad << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2BeforeModElem[i].m_Ecc << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2BeforeModElem[i].m_I* AsCRadToDeg << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2BeforeModElem[i].m_RAAN * AsCRadToDeg << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2BeforeModElem[i].m_ArgPeri* AsCRadToDeg << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2BeforeModElem[i].m_TrueA* AsCRadToDeg << endl << endl
				<< "MJ2000�����������:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2ModElem[i].m_PeriRad << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2ModElem[i].m_Ecc << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2ModElem[i].m_I *AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2ModElem[i].m_RAAN *AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2ModElem[i].m_ArgPeri *AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2ModElem[i].m_TrueA *AsCRadToDeg << endl << endl
				<< "----------" << endl
				<< "�����α��" << endl
				<< "----------" << endl << endl
				<< "�����Ԫ��UTCG��:   " << endl << endl
				<< setw(18) << setprecision(15) << timestr3 << endl << endl
				<< "������� MJ2000��m/s��:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3DeltaV[i][0] << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3DeltaV[i][1] << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3DeltaV[i][2] << endl << endl
				<< "ת��ǰMJ2000�����������:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3BeforeModElem[i].m_PeriRad << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3BeforeModElem[i].m_Ecc << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3BeforeModElem[i].m_I* AsCRadToDeg << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3BeforeModElem[i].m_RAAN * AsCRadToDeg << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3BeforeModElem[i].m_ArgPeri* AsCRadToDeg << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3BeforeModElem[i].m_TrueA* AsCRadToDeg << endl << endl
				<< "MJ2000�����������:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3ModElem[i].m_PeriRad << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3ModElem[i].m_Ecc << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3ModElem[i].m_I *AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3ModElem[i].m_RAAN *AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3ModElem[i].m_ArgPeri *AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3ModElem[i].m_TrueA *AsCRadToDeg << endl << endl
				<< "------------" << endl
				<< "�����״̬" << endl
				<< "------------" << endl << endl
				<< "������Ԫ��UTCG��:   " << endl << endl
				<< setw(18) << setprecision(15) << timeend << endl << endl
				<< "MJ2000�����������:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LDOModElem[i].m_PeriRad << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LDOModElem[i].m_Ecc << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LDOModElem[i].m_I*AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LDOModElem[i].m_RAAN*AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LDOModElem[i].m_ArgPeri*AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LDOModElem[i].m_TrueA*AsCRadToDeg << endl << endl
				<< "���ٶ�������m/s��:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TotalImpulse[i] << endl << endl
				<< "����deg��:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.epsi[i] * AsCRadToDeg << endl << endl
				<< "���߼нǣ�deg��:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.perierror[i] * AsCRadToDeg << endl << endl;
		}
	}
	outfile.close();
	return true;
}


















