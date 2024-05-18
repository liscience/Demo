#include "FreeReturnOrbit.h"
#include <Eigen/SVD>
#include <Eigen/Core>
#include<iostream>
#include<fstream>
#include<cmath>
#include<iomanip>
using namespace std;

//构造函数
CFreeROInit::CFreeROInit()
{
	;
}
//析构函数
CFreeROInit::~CFreeROInit()
{
	;
}



//**********************************************************
///NR方法的自由返回轨道初步计算主流程函数
///@Author   Li Zeyue
///@Date	 2022.7.30
///@Input    
///			 环月目标轨道修正轨道根数 LLOModOrbElem
///          设计变量 x
///          右函数目标计算结果 Y_tar
///          误差限 Out_Ylim
///          设计变量参考区间 int_interval
///          设计变量迭代范围 Iter_interval
///          区间比例因子 k
///          迭代次数 Iternum
///@Output
///          是否收敛的布尔值 IsCon
//***********************************************************
void CFreeROInit::NRFreeRTOInit(const CModOrbElem& LLOModOrbElem, CVector<double>& x,
	CVector<double> Y_tar, CVector<double> Out_Ylim, CVector<double> int_interval,
	CMatrix<double> Iter_interval, CVector<double> k, int Iternum, bool&IsCon)
{

	double jd1, jd2;                 //儒略日
	double deltaV;                   //速度增量
	CModOrbElem TLI_LTOModElem;      //出发时刻转移轨道根数
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
	//搜索开始和结束时间
	AsTimeToJD(mPRLTime[0], jd1);
	AsTimeToJD(mPRLTime[1], jd2);
	int Optxnum = x.GetSize();
	mLLOModElem = LLOModOrbElem;
	for (double i = jd1; i <= jd2; i = i + mStep*AsCSecToDay)
	{
		AsJDToTime(i, mLOITime);
		mLLOModElemNew = mLLOModElem;
		//生成新的轨道根数
		AsModOrbElemToCart(mLLOModElem, AsCMoonGrav, pos, vel);
		propSelenocentricElem.m_Epoch = i;
		prop.LunarOrbitPropVarStep(propSelenocentricElem, mStep, 1e-13, 1, pos, vel, timeArray, posArray, velArray);
		AsCartToModOrbElem(pos, vel, AsCMoonGrav, mLLOModElem);
		NRMethod(&CFreeROInit::ConJud, &CFreeROInit::SenCalFun, &CFreeROInit::SenCal,
			&CFreeROInit::SenCalFun, x, Y_tar, Out_Ylim,
			int_interval, Iter_interval, k, Iternum, IsCon);
		if (IsCon != 0)
		{
			//导出结果
			mFreeReturnOrbitInitCalOutput.PRLTime.push_back(mLOITime);
			double deltat;
			deltat = AsTrueToTPP(mLTOTrueA, mLSIModElemEJ2000.m_PeriRad / (1 - mLSIModElemEJ2000.m_Ecc), mLSIModElemEJ2000.m_Ecc, AsCEarthGrav);

			mTLITime.AddSec(deltat);
			mFreeReturnOrbitInitCalOutput.TLILaunchTime.push_back(mTLITime);

			//记录进入奔月轨道倾角
			mFreeReturnOrbitInitCalOutput.LTOInc.push_back(mLSIModElemEJ2000.m_I);
			//记录伪参数
			mFreeReturnOrbitInitCalOutput.PRLPseudoLat.push_back(mPRLPseudoLat);
			mFreeReturnOrbitInitCalOutput.PRLPseudoLong.push_back(mPRLPseudoLon);
			mFreeReturnOrbitInitCalOutput.PRLEccentric.push_back(mPRLEccentric);
			mFreeReturnOrbitInitCalOutput.PRLVelInc.push_back(mPRLVelInc);
			//记录月心轨道参数
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
			//计算真空近地点处的位置速度
			CalHPLOI2VCP(mLOITime, LCO_PRLPosMJ2000, LCO_PRLVelMJ2000, tempTime, PosEJ2000, VelEJ2000, timeArray, posArray, velArray);

			//计算ROI的地理经纬度
			CalROIECFLatLon(timeArray, posArray, velArray, mLOITime, mAtmospH, timeLOI2ROI, mROIECFLat, mROIECFLon);
			tempTime = mLOITime;
			tempTime.AddSec(timeLOI2ROI);
			mFreeReturnOrbitInitCalOutput.ROITime.push_back(tempTime);
			mFreeReturnOrbitInitCalOutput.ROILat.push_back(mROIECFLat);
			mFreeReturnOrbitInitCalOutput.ROILon.push_back(mROIECFLon);

			//奔月地心段
			if (fabs(AsNegPIToPI(mLSIModElemEJ2000.m_ArgPeri)) <= AsCHalfPI)
			{
				mFreeReturnOrbitInitCalOutput.IsLTODesend.push_back(false);//奔月地心段升轨
			}
			else
			{
				mFreeReturnOrbitInitCalOutput.IsLTODesend.push_back(true);//奔月地心段降轨
			}
			//奔月月心段
			if (fabs(AsNegPIToPI(mLCO_PRLModElem.m_ArgPeri)) <= AsCHalfPI)
			{
				mFreeReturnOrbitInitCalOutput.IsLDODesend.push_back(false);//月心段升轨
			}
			else
			{
				mFreeReturnOrbitInitCalOutput.IsLDODesend.push_back(true);//月心段降轨
			}
			//返回地心段
			if (fabs(AsNegPIToPI(mLSOModElemEJ2000.m_ArgPeri)) <= AsCHalfPI)
			{
				mFreeReturnOrbitInitCalOutput.IsRTODesend.push_back(false);//地心段升轨
			}
			else
			{
				mFreeReturnOrbitInitCalOutput.IsRTODesend.push_back(true);//地心段降轨
			}
			goto OUTLOOP;
		}
	OUTLOOP:
		continue;
	}
}

//**********************************************************
///NR方法的自由返回轨道高精度计算主流程函数
///@Author   Li Zeyue
///@Date	 2022.7.30
///@Input    
///			 环月目标轨道修正轨道根数 LLOModOrbElem
///          设计变量 x
///          右函数目标计算结果 Y_tar
///          误差限 Out_Ylim
///          设计变量参考区间 int_interval
///          设计变量迭代范围 Iter_interval
///          区间比例因子 k
///          迭代次数 Iternum
///@Output
///          是否收敛的布尔值 IsCon
//***********************************************************
void CFreeROInit::NRFreeRTOHP(const CModOrbElem& LLOModOrbElem, CVector<double>& x,
	CVector<double> Y_tar, CVector<double> Out_Ylim, CVector<double> int_interval,
	CMatrix<double> Iter_interval, CVector<double> k, int Iternum, bool&IsCon)
{

	double jd1, jd2;                 //儒略日
									 //搜索开始和结束时间
	AsTimeToJD(mPRLTime[0], jd1);
	AsTimeToJD(mPRLTime[1], jd2);
	int size;
	size = (jd2 - jd1) / (mStep*AsCSecToDay) + 1;
	//设置高精度输出结构体大小

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
	double       deltV;                 //制动速度
	double       reentryGamma;          //再入角

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
	//星历数据
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
		//生成新的轨道根数
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
			//记录优化结果 
			Cal1Impulse(mTLIModOrbElem, 0, deltV);
			tempTime = mLOITime;
			double deltat;
			//deltat = AsTrueToTPP(mLTOTrueA, mTLIModOrbElem.m_PeriRad / (1 - mTLIModOrbElem.m_Ecc), mTLIModOrbElem.m_Ecc, AsCEarthGrav);
			AsModOrbElemToCart(mLCO_PRLModElem, AsCMoonGrav, PRL_LCOPosMJ2000, PRL_LCOVelMJ2000);
			// 高精度模型由近月点伪状态逆向求解近地出发时刻状态
			CalHPLOI2TLI(mLOITime, PRL_LCOPosMJ2000, PRL_LCOVelMJ2000, tempTime, PosEJ2000, VelEJ2000, timeArray, posArray, velArray);
			tempTime = mLOITime;

			//申请动态内存			
			CModOrbElem	 tempTLI_LTOModeElem;

			int number = timeArray.size() / 2;
			int counter = 0;
			double *p_time = new double[number];
			double *p_trueanal = new double[number];

			//求解TLI时刻
			for (unsigned int i = timeArray.size() - number; i < timeArray.size(); i++)
			{
				//各时刻位置
				p_time[counter] = timeArray[i];
				AsCartToModOrbElem(posArray[i], velArray[i], AsCEarthGrav, tempTLI_LTOModeElem);
				p_trueanal[counter] = AsNegPIToPI(tempTLI_LTOModeElem.m_TrueA);
				counter++;
			}
			//插值获得到达载入点的时间
			AsInterpPol1D(number, p_trueanal, p_time, 7, false, mLTOTrueA, deltat);
			AsInterpPosVel(timeArray, posArray, velArray, 7, false, deltat, PosEJ2000, VelEJ2000);
			tempTime.AddSec(deltat);
			// TLI时刻 (UTC) 
			mFreeRTONRCalOutput.TLILaunchTime[num] = tempTime;
			// TLI时刻，LTO的轨道修正根数

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
			//非全相位
			//mperierror = (asin(sin(mLCO_PRLModElem.m_RAAN - RAAN)*sin(LLOInc) / sin(mepsi)) - mLCO_PRLModElem.m_ArgPeri);

			mFreeRTONRCalOutput.epsi[num] = mepsi;
			mFreeRTONRCalOutput.perierror[num] = asin(mperierror);
			mFreeRTONRCalOutput.ReentryGamma[num] = reentryGamma;

			AsModOrbElemToCart(mLCO_PRLModElem, AsCMoonGrav, PRL_LCOPosMJ2000, PRL_LCOVelMJ2000);
			//计算VCP处的位置、速度
			CalHPLOI2VCP(mLOITime, PRL_LCOPosMJ2000, PRL_LCOVelMJ2000, tempTime1, PosEJ2000, VelEJ2000, timeArray, posArray, velArray);
			if (timeArray.size() == 0)
			{
				//超出误差限
				cout << "时间序列为0" << endl;
			}
			//计算ROI的地理经纬度
			CalROIECFLatLon(timeArray, posArray, velArray, mLOITime, mAtmospH, timeLOI2ROI, mROIECFLat, mROIECFLon);
			tempTime = mLOITime;
			tempTime.AddSec(timeLOI2ROI);
			//计算真空近地点和返回参数
			mFreeRTONRCalOutput.ROITime[num] = tempTime;
			mFreeRTONRCalOutput.ROIECFLat[num] = mROIECFLat;
			mFreeRTONRCalOutput.ROIECFLon[num] = mROIECFLon;

			CCoord3 ROIPos, ReturnPos;
			CLatLonRad ROIECF(mROIECFLat, mROIECFLon, AsCEarthRadius);
			double    angROI2Return;

			//计算再入轨道上返回场纬度对应的地理经度
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

			//计算升交点经度
			CCoord TLIpos, TLIvel;
			CCoord TLIposECF, TLIvelECF;
			double jdp1, jdp2;
			CModOrbElem TLI_LTOModElem;
			CMatrix<double> J2000ToECF(3, 3);
			TLI_LTOModElem = mTLIModOrbElem;
			//地固系下的根数
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

			////计算升交点赤经
			double siteLon, TLIraan;
			//确定TLI时刻发射场共面轨道的升交点经度、象限
			DeterminPhase(mFreeRTONRCalOutput.TLILaunchTime[num], mTLIModOrbElem, findBackwardtimeFunc.mLaunchSite,
				findBackwardtimeFunc.mIsDescend, findBackwardtimeFunc.mflag, siteLon);
			//初始化时间、根数
			findBackwardtimeFunc.mTLITime = mFreeRTONRCalOutput.TLILaunchTime[num];
			findBackwardtimeFunc.mLTIOrbElem = mTLIModOrbElem;
			TLIraan = mTLIModOrbElem.m_RAAN;
			//发射升交点在1、2象限、TLI升交点在3、4象限
			if (findBackwardtimeFunc.mflag <3)
			{
				TLIraan = AsNegPIToPI(TLIraan);
				siteLon = AsNegPIToPI(siteLon);
			}
			if (siteLon - TLIraan < AsCPI && siteLon - TLIraan > 0)
			{
				//求解回退时间,解范围为0~12小时
				AsNLEBisection(findBackwardtimeFunc, -43200, 0, 1e-6, backwardTime);
			}
			else
			{
				//考虑迭代方向，将象限反置
				if (findBackwardtimeFunc.mflag <= 2)
				{
					findBackwardtimeFunc.mflag += 2;
				}
				else
				{
					findBackwardtimeFunc.mflag -= 2;
				}
				//求解回退时间,解范围为12~24小时
				AsNLEBisection(findBackwardtimeFunc, -86400, -43200, 1e-6, backwardTime);
			}

			//求解共面时刻的相位差
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

// 利用Eigen库，采用SVD分解的方法求解矩阵伪逆，默认误差er为0
Eigen::MatrixXd pinv_eigen_based(Eigen::MatrixXd& origin, const float er = 0) {
	// 进行svd分解
	Eigen::JacobiSVD<Eigen::MatrixXd> svd_holder(origin,
		Eigen::ComputeThinU |
		Eigen::ComputeThinV);
	// 构建SVD分解结果
	Eigen::MatrixXd U = svd_holder.matrixU();
	Eigen::MatrixXd V = svd_holder.matrixV();
	Eigen::MatrixXd D = svd_holder.singularValues();

	// 构建S矩阵
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
///NR方法主流程函数（支持初步计算和高精度计算的混合操作）
///@Author   Li Zeyue
///@Date	 2022.7.30
///@Input    
///			 收敛性判断函数指针 Confunc类成员
///			 右函数计算函数指针 Calfunc类成员
///			 敏感度矩阵计算函数指针 SenMa类成员
///          设计变量 x
///          右函数目标计算结果 Y_tar
///          误差限 Out_Ylim
///          设计变量参考区间 int_interval
///          设计变量迭代范围 Iter_interval
///          区间比例因子 k
///          迭代次数 Iternum
///@Output
///          是否收敛的布尔值 IsCon
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
		Optx[i] = x[i];//将迭代初值赋值给Optx

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
			//SenInv = sen_result.Inv();//得到当前点的数值导数逆矩阵
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
		NRStep = SenInv.operator*(yback);//逆矩阵乘以迭代参数
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
///NR方法的自由返回轨道计算收敛性判断
///@Author   Li Zeyue
///@Date	 2022.7.30
///@Input    
///			 右函数实际计算结果 Out_Y
///          右函数目标计算结果 Y_tar
///          误差限 Out_Ylim
///@Output
///          是否收敛的布尔值 IsCon
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
///基于近月点伪参数的右函数计算函数（双二体初步计算）
///@Author   Li Zeyue
///@Date	 2022.7.30
///@Input    
///			 设计变量 x
///@Output
///          右函数计算结果 result
//***********************************************************
void CFreeROInit::SenCalFun(CVector<double>& x, CVector<double>& result)
{
	//计算目标值
	double timeofFlight;           //飞行时间
	double tLSI2LOI;               //LSI到LOI的转移时间
	double timeLSO2VCP;            //LSO到VCP的转移时间
	double timeLSO2ROI;            //LSO到ROI的转移时间
	double rr, DELTA;              //计算真空近地点的参数
								   //CVector<double> observations;
								   //observations.Resize(12);
	CCoord3 LCO_PRLPosMJ2000, LCO_PRLVelMJ2000;//近月点处LCO轨道的位置速度
	CCTime LOITime;                //近月加速时刻

								   //将优化的参数传给计算优化目标的功能函数
	mPRLPseudoLat = x[0];
	mPRLPseudoLon = x[1];
	mPRLEccentric = x[2];
	mPRLVelInc = x[3];
	//mPRLPseudoLat = -0.016996450218444769;
	//mPRLPseudoLon = 3.1173658551551853;
	//mPRLEccentric = 1.5113342651185437;
	//mPRLVelInc =	0.15867919986113038;
	LOITime = mLOITime;
	//由近月点伪状态逆向求解近地出发时刻状态（双二体模型）
	CalInitPseudoPRL2LSI(LOITime, mLCO_PRLHeight, mPRLPseudoLat, mPRLPseudoLon, mPRLEccentric, mPRLVelInc, mLSIModElemEJ2000);
	// 计算从月球影响球到近地点的时间
	CalTimeFromTrueA(mLSIModElemEJ2000, 0, AsCEarthGrav, timeofFlight);
	timeofFlight = fabs(timeofFlight);
	PerilunePseudoState2MJ2000(LOITime, mLCO_PRLHeight, mPRLPseudoLat, mPRLPseudoLon, mPRLEccentric, mPRLVelInc,
		LCO_PRLPosMJ2000, LCO_PRLVelMJ2000);
	//计算近月点制动前的轨道根数
	AsCartToModOrbElem(LCO_PRLPosMJ2000, LCO_PRLVelMJ2000, AsCMoonGrav, mLCO_PRLModElem);
	CalTimeLSI2LOI(mLCO_PRLModElem, tLSI2LOI);
	mTofTemp = timeofFlight + tLSI2LOI;
	//出发时刻
	mTLITime = mLOITime;
	mTLITime.AddSec(-mTofTemp);
	//由近月点伪状态正向求解返回出发时刻状态（双二体模型）
	CalInitPseudoPRL2LSO(LOITime, mLCO_PRLHeight, mPRLPseudoLat, mPRLPseudoLon,
		mPRLEccentric, mPRLVelInc, mLSOModElemEJ2000);
	//LSO到VCP时间
	CalTimeFromTrueA(mLSOModElemEJ2000, AsCTwoPI - 0.001, AsCEarthGrav, timeLSO2VCP);
	timeLSO2VCP = fabs(timeLSO2VCP);

	//LSO到ROI时间
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

	//求解真空近地点高度，一般为5~51km
	rr = mAtmospH + AsCEarthRadius;
	DELTA = AsSqr(rr)*(1 - (1 - AsSqr(mLSOModElemEJ2000.m_Ecc))*(AsSqr(tan(mReentryGamma)) + 1));

	if (mLSOModElemEJ2000.m_Ecc < 0.96)
	{
		//未被捕获
		mVCPH = 50000;
	}
	else
	{
		//求解真空近地点高度
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
///基于近月点轨道根数的右函数计算函数（双二体初步计算）
///@Author   Li Zeyue
///@Date	 2022.7.30
///@Input    
///			 设计变量 x
///@Output
///          右函数计算结果 result
//***********************************************************
void CFreeROInit::SenCalFunElements(CVector<double>& x, CVector<double>& result)
{
	//计算目标值

	double timeofFlight;           //飞行时间
	double tLSI2LOI;               //LSI到LOI的转移时间
	double timeLSO2VCP;            //LSO到VCP的转移时间
	double timeLSO2ROI;            //LSO到ROI的转移时间
	double rr, DELTA;              //计算真空近地点的参数
								   //CVector<double> observations;
								   //observations.Resize(12);
	CCoord3 LCO_PRLPosMJ2000, LCO_PRLVelMJ2000;//近月点处LCO轨道的位置速度
	CCTime LOITime;                //近月加速时刻

								   //将优化的参数传给计算优化目标的功能函数
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
	//由近月点月心段轨道根数逆向求解近地出发时刻状态（双二体模型）
	CalInitElementsPRL2LSI(LOITime, mLCO_PRLHeight, PRLPeriluneEccentric, PRLPeriluneInc,
		PRLPeriluneRAAN, PRLPeriluneLAOP, mLSIModElemEJ2000);
	// 计算从月球影响球到近地点的时间
	CalTimeFromTrueA(mLSIModElemEJ2000, 0, AsCEarthGrav, timeofFlight);
	timeofFlight = fabs(timeofFlight);
	//计算近月点制动前的轨道根数
	CalTimeLSI2LOI(mLCO_PRLModElem, tLSI2LOI);
	mTofTemp = timeofFlight + tLSI2LOI;
	//出发时刻
	mTLITime = mLOITime;
	mTLITime.AddSec(-mTofTemp);
	//由月心段轨道根数正向求解返回出发时刻状态（双二体模型）
	CalInitElementsPRL2LSO(LOITime, mLCO_PRLHeight, PRLPeriluneEccentric, PRLPeriluneInc,
		PRLPeriluneRAAN, PRLPeriluneLAOP, mLSOModElemEJ2000);
	//LSO到VCP时间
	CalTimeFromTrueA(mLSOModElemEJ2000, AsCTwoPI - 0.001, AsCEarthGrav, timeLSO2VCP);
	timeLSO2VCP = fabs(timeLSO2VCP);

	//LSO到ROI时间
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

	//求解真空近地点高度，一般为5~51km
	rr = mAtmospH + AsCEarthRadius;
	DELTA = AsSqr(rr)*(1 - (1 - AsSqr(mLSOModElemEJ2000.m_Ecc))*(AsSqr(tan(mReentryGamma)) + 1));

	if (mLSOModElemEJ2000.m_Ecc < 0.96)
	{
		//未被捕获
		mVCPH = 50000;
	}
	else
	{
		//求解真空近地点高度
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
///基于近月点伪参数的右函数计算函数（高精度计算）
///@Author   Li Zeyue
///@Date	 2022.7.30
///@Input    
///			 设计变量 x
///@Output
///          右函数计算结果 result
//***********************************************************
void CFreeROInit::SenHPCalFun(CVector<double>& x, CVector<double>& result)
{
	//伪参数
	double PRLPseudoLatNew;
	double PRLPseudoLongNew;
	double PRLPeriluneVelIncNew;
	double PRLPeriluneEccentricNew;
	double rr, DELTA;              //计算真空近地点的参数
								   //double epsi, perierror;
								   //double reentryGamma;                  //再入角
	CCoord3 PRL_LCOPosMJ2000;             //LCO轨道近月点位置
	CCoord3 PRL_LCOVelMJ2000;             //LCO轨道近月点速度
	CCoord3 PosEJ2000, VelEJ2000;
	CCTime LOITime;                       //近月点机动时间
	CCTime tempTime;
	//星历数据
	std::vector<double> timeArray;
	std::vector<CCoord> posArray;
	std::vector<CCoord> velArray;

	/////////////将优化的参数传给计算优化目标的功能函数////////////
	PRLPseudoLatNew = x[0];
	PRLPseudoLongNew = x[1];
	PRLPeriluneEccentricNew = x[2];
	PRLPeriluneVelIncNew = x[3];
	LOITime = mLOITime;

	// 由伪参数获得位置速度
	PerilunePseudoState2MJ2000(LOITime, mLCO_PRLHeight, PRLPseudoLatNew, PRLPseudoLongNew, PRLPeriluneEccentricNew,
		PRLPeriluneVelIncNew, PRL_LCOPosMJ2000, PRL_LCOVelMJ2000);
	AsCartToModOrbElem(PRL_LCOPosMJ2000, PRL_LCOVelMJ2000, AsCMoonGrav, mLCO_PRLModElem);
	// 高精度模型由近月点伪状态逆向求解近地出发时刻状态
	CalHPLOI2TLI(LOITime, PRL_LCOPosMJ2000, PRL_LCOVelMJ2000, tempTime, PosEJ2000, VelEJ2000, timeArray, posArray, velArray);
	AsCartToModOrbElem(PosEJ2000, VelEJ2000, AsCEarthGrav, mTLIModOrbElem);
	mTof = (tempTime - LOITime).GetTotalSec();
	// 高精度模型由近月点正向求解返回近地时刻状态
	// 精确求解某条自由返回轨道时采用下面的函数，如果发现轨道倾角值不是目标值，需要调整近月点时间
	CalHPLOI2VCP(LOITime, PRL_LCOPosMJ2000, PRL_LCOVelMJ2000, tempTime, PosEJ2000, VelEJ2000, timeArray, posArray, velArray);
	AsCartToModOrbElem(PosEJ2000, VelEJ2000, AsCEarthGrav, mVCPModOrbElem);
	//求解真空近地点高度，一般为5~51km
	rr = mAtmospH + AsCEarthRadius;
	DELTA = AsSqr(rr)*(1 - (1 - AsSqr(mLSOModElemEJ2000.m_Ecc))*(AsSqr(tan(mReentryGamma)) + 1));

	if (mLSOModElemEJ2000.m_Ecc < 0.96)
	{
		//未被捕获
		mVCPH = 50000;
	}
	else
	{
		//求解真空近地点高度
		mVCPH = (rr + sqrt(DELTA)) / (1 + mLSOModElemEJ2000.m_Ecc) / (1 + AsSqr(tan(mReentryGamma))) - AsCEarthRadius;
	}

	//根据经验固定加权函数
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
///基于近月点轨道根数的右函数计算函数（高精度计算）
///@Author   Li Zeyue
///@Date	 2022.9.7
///@Input    
///			 设计变量 x
///@Output
///          右函数计算结果 result
//***********************************************************
void CFreeROInit::SenHPCalFunElement(CVector<double>& x, CVector<double>& result)
{
	//伪参数
	double PRLPeriluneEccentric;
	double PRLPeriluneInc;
	double PRLPeriluneRAAN;
	double PRLPeriluneLAOP;

	double rr, DELTA;              //计算真空近地点的参数
	//double epsi, perierror;
	//double reentryGamma;                  //再入角
	CCoord3 PRL_LCOPosMJ2000;             //LCO轨道近月点位置
	CCoord3 PRL_LCOVelMJ2000;             //LCO轨道近月点速度
	CCoord3 PosEJ2000, VelEJ2000;
	CCTime LOITime;                       //近月点机动时间
	CCTime tempTime;
	//星历数据
	std::vector<double> timeArray;
	std::vector<CCoord> posArray;
	std::vector<CCoord> velArray;

	/////////////将优化的参数传给计算优化目标的功能函数////////////
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
	// 高精度模型由近月点伪状态逆向求解近地出发时刻状态
	CalHPLOI2TLI(LOITime, PRL_LCOPosMJ2000, PRL_LCOVelMJ2000, tempTime, PosEJ2000, VelEJ2000, timeArray, posArray, velArray);
	AsCartToModOrbElem(PosEJ2000, VelEJ2000, AsCEarthGrav, mTLIModOrbElem);
	mTof = (tempTime - LOITime).GetTotalSec();
	// 高精度模型由近月点正向求解返回近地时刻状态
	// 精确求解某条自由返回轨道时采用下面的函数，如果发现轨道倾角值不是目标值，需要调整近月点时间
	CalHPLOI2VCP(LOITime, PRL_LCOPosMJ2000, PRL_LCOVelMJ2000, tempTime, PosEJ2000, VelEJ2000, timeArray, posArray, velArray);
	AsCartToModOrbElem(PosEJ2000, VelEJ2000, AsCEarthGrav, mVCPModOrbElem);
	//求解真空近地点高度，一般为5~51km
	rr = mAtmospH + AsCEarthRadius;
	DELTA = AsSqr(rr) * (1 - (1 - AsSqr(mLSOModElemEJ2000.m_Ecc)) * (AsSqr(tan(mReentryGamma)) + 1));

	if (mLSOModElemEJ2000.m_Ecc < 0.96)
	{
		//未被捕获
		mVCPH = 50000;
	}
	else
	{
		//求解真空近地点高度
		mVCPH = (rr + sqrt(DELTA)) / (1 + mLSOModElemEJ2000.m_Ecc) / (1 + AsSqr(tan(mReentryGamma))) - AsCEarthRadius;
	}

	//根据经验固定加权函数
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
///自由返回轨道计算参数敏感度矩阵（前向差分方法）
///@Author   Li Zeyue
///@Date	 2022.7.30
///@Input    
///			 右函数计算函数指针 Calfunc类成员
///          设计变量 x
///          设计变量参考区间 int_interval
///          区间比例因子 k
///@Output
///          参数敏感性矩阵 sen_result
///          右函数计算结果 result
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
/// 近月点伪状态转化到月心J2000位置速度
/// @Author	He Boyong   Zhou Wanmeng
/// @Date	2015-01-21  2018.5.12
/// @Input
/// @Param	timePRL						近月点时刻，UTCG
/// @Param	PRLHeight					近月点高度,m
/// @Param	PRLPerilunePseudoLat		近月点在月心瞬时惯性系中方位角纬度，rad
/// @Param	PRLPerilunePseudoLong		近月点在月心瞬时惯性系中方位角经度，rad
/// @Param	PRLPeriluneEccentric  		近月点LCO轨道偏心率，1.3，等价与近月点速度（m/s）一般在2.55e3 m/s
/// @Param	PRLPeriluneVelInc			近月点伪倾角，可以在pi±0.20943951，rad = 180±12，deg
/// @Output
/// @Param	posMJ2000				    月心J2000位置，m
/// @Param	velMJ2000				    月心J2000速度，m/s
//***********************************************************
void CFreeROInit::PerilunePseudoState2MJ2000(const CCTime& timePRL, const double& PRLHeight, const double& PRLPerilunePseudoLat,
	const double& PRLPerilunePseudoLong, const double& PRLPeriluneEccentric, const double& PRLPeriluneVelInc,
	CCoord3& posMJ2000, CCoord3& velMJ2000)
{
	CCoord3 posMPI, velMPI;				//月心近月点坐标系
	CCoord3 posLVLH, velLVLH;			//月心轨道瞬时惯性系
	CMatrix<double> LVLHtoMPI(3, 3);    //LVLH到MPI坐标系
	CMatrix<double> LVLHtoMJ2000(3, 3); //LVLH到MJ2000坐标系
	double semiMajor, PRLPeriluneVel;   //半长轴、近月点速度
	double prlJd;

	//半长轴
	semiMajor = (PRLHeight + AsCMoonRadius) / (1 - PRLPeriluneEccentric);
	//近月点速度
	PRLPeriluneVel = sqrt(AsCMoonGrav*(2.0 / (PRLHeight + AsCMoonRadius) - 1.0 / semiMajor));
	//近月点位置、速度
	posMPI[0] = PRLHeight + AsCMoonRadius;
	posMPI[1] = 0;
	posMPI[2] = 0;
	velMPI[0] = 0;
	velMPI[1] = PRLPeriluneVel*cos(PRLPeriluneVelInc);
	velMPI[2] = PRLPeriluneVel*sin(PRLPeriluneVelInc);

	//轨道瞬时惯性系MOI到近月点坐标系MPI的转化矩阵
	AsEulerToMtx(CEuler(PRLPerilunePseudoLong, -PRLPerilunePseudoLat, 0), 321, LVLHtoMPI);
	posLVLH = LVLHtoMPI.Transpose()*posMPI;
	velLVLH = LVLHtoMPI.Transpose()*velMPI;
	//月地转移出发时刻对应儒略日
	AsTimeToJD(timePRL, prlJd);
	//生成近月点坐标系到轨道瞬时惯性系的转化矩阵
	MoonLVLHToJ2000(timePRL, LVLHtoMJ2000);

	//转化到月心J2000坐标系
	posMJ2000 = LVLHtoMJ2000*posLVLH;
	velMJ2000 = LVLHtoMJ2000*velLVLH;
}
//***********************************************************
/// 由近月点伪状态逆向求解近地出发时刻状态（双二体模型）
/// @Author	He Boyong   Zhou Wanmeng
/// @Date	2015-01-21  2018.5.12
/// @Input
/// @Param	timePRL						近月点时刻，UTCG
/// @Param	PRLHeight					近月点高度,m
/// @Param	PRLPerilunePseudoLat		近月点在月心瞬时惯性系中方位角纬度，rad
/// @Param	PRLPerilunePseudoLong		近月点在月心瞬时惯性系中方位角经度，rad
/// @Param	PRLPeriluneEccentric  		近月点LCO轨道偏心率，1.3，等价与近月点速度（m/s）一般在2.55e3 m/s
/// @Param	PRLPeriluneVelInc			近月点伪倾角，可以在pi±0.20943951，rad = 180±12，deg
/// @Output
/// @Param	LSIModElemEJ2000      进入月球影响球时刻EJ2000的修正轨道根数
/// @Return		
//***********************************************************
void CFreeROInit::CalInitPseudoPRL2LSI(const CCTime& timePRL, const double& PRLHeight, const double& PRLPerilunePseudoLat, const double& PRLPerilunePseudoLong,
	const double& PRLPeriluneEccentric, double PRLPeriluneVelInc, CModOrbElem& LSIModElemEJ2000)
{
	CCoord3 posMJ2000, velMJ2000;
	CModOrbElem PRL_LCOModOrbElem;

	//由伪状态转换为近月点MJ2000位置速度
	PerilunePseudoState2MJ2000(timePRL, PRLHeight, PRLPerilunePseudoLat,
		PRLPerilunePseudoLong, PRLPeriluneEccentric, PRLPeriluneVelInc, posMJ2000, velMJ2000);
	//由位置速度计算近月点LCO修正轨道根数
	AsCartToModOrbElem(posMJ2000, velMJ2000, AsCMoonGrav, PRL_LCOModOrbElem);
	//由PRL_LCOModOrbElem计算月球影响球入口点的EJ2000轨道参数
	CalInitMJ2000PRL2LSI(timePRL, PRL_LCOModOrbElem, PRL_LCOModOrbElem.m_ArgPeri + PRL_LCOModOrbElem.m_TrueA,
		PRL_LCOModOrbElem.m_Ecc, LSIModElemEJ2000);
}
void CFreeROInit::CalInitElementsPRL2LSI(const CCTime& timePRL, const double& PRLHeight, const double& PRLPeriluneEccentric, const double& PRLPeriluneInc,
	const double& PRLPeriluneRAAN, double PRLPeriluneLAOP, CModOrbElem& LSIModElemEJ2000)
{
	CCoord3 posMJ2000, velMJ2000;
	CModOrbElem PRL_LCOModOrbElem;
	//近月点LCO修正轨道根数
	PRL_LCOModOrbElem.m_PeriRad = PRLHeight + AsCMoonRadius;
	PRL_LCOModOrbElem.m_Ecc = PRLPeriluneEccentric;
	PRL_LCOModOrbElem.m_I = PRLPeriluneInc;
	PRL_LCOModOrbElem.m_RAAN = PRLPeriluneRAAN;
	PRL_LCOModOrbElem.m_ArgPeri = PRLPeriluneLAOP;
	PRL_LCOModOrbElem.m_TrueA = 0;
	//由PRL_LCOModOrbElem计算月球影响球入口点的EJ2000轨道参数
	CalInitMJ2000PRL2LSI(timePRL, PRL_LCOModOrbElem, PRL_LCOModOrbElem.m_ArgPeri + PRL_LCOModOrbElem.m_TrueA,
		PRL_LCOModOrbElem.m_Ecc, LSIModElemEJ2000);
}
/// 高精度模型由近月点正向求解返回近地时刻状态
/// @Author	He Boyong   Zhou Wanmeng
/// @Date	2015-04-12  2018.7.3
/// @Input
/// @Param	timePRL		            近月点时刻,UTCG
/// @Param	PRL_LCOPosMJ2000		月心J2000坐标系双曲线轨道近月点位置,m
/// @Param	PRL_LCOVelMJ2000		月心J2000坐标系双曲线轨道近月点速度,m/s
/// @Output
/// @Param	time		            近地点时刻,UTCG
/// @Param	PosEJ2000				地心J2000坐标系LTO轨道近地点位置,m
/// @Param	VelJ2000				地心J2000坐标系LTO轨道近地点速度,m/s
/// @Param  timeArray               历史数据
/// @Param  posArray
/// @Param  velArray
/// @Return	
//***********************************************************
bool CFreeROInit::CalHPLOI2VCP(const CCTime& timePRL, const CCoord3& PRL_LCOPosMJ2000, const CCoord3& PRL_LCOVelMJ2000,
	CCTime& time, CCoord3& PosEJ2000, CCoord3& VelEJ2000, std::vector<double>& timeArray,
	std::vector<CCoord3>& posArray, std::vector<CCoord3>& velArray)
{
	double PRLJd;                                //PRL时刻
	CCTime TT;
	CModOrbElem PRL_LCOModOrbElem;               //近月点参数的修正轨道根数
	CCoord3 lunarPosEJ2000, lunarVelEJ2000;      //月球在地心J2000坐标系的参数
	CCoord3 PRL_LCOPosEJ2000, PRL_LCOVelEJ2000;  //近月点相对地球的EJ2000位置速度
	bool    isSuccess;

	//数据清空
	timeArray.resize(0);
	posArray.resize(0);
	velArray.resize(0);

	//计算PRL时刻
	AsUTCToTT(timePRL, TT);
	AsTimeToJD(TT, PRLJd);

	//月球在地心J2000坐标系的参数
	AsJplPlanetaryEph(PRLJd, AsEJplMoon, AsEJplEarth, lunarPosEJ2000, lunarVelEJ2000);
	//近月点相对地球的EJ2000位置速度
	PRL_LCOPosEJ2000 = lunarPosEJ2000 + PRL_LCOPosMJ2000;
	PRL_LCOVelEJ2000 = lunarVelEJ2000 + PRL_LCOVelMJ2000;
	//LCO的修正轨道根数EJ2000
	AsCartToModOrbElem(PRL_LCOPosEJ2000, PRL_LCOVelEJ2000, AsCEarthGrav, PRL_LCOModOrbElem);

	///////////////配置积分停止条件////////////////
	CStateTrigger stopperigee;
	stopperigee.m_TripType = 3;
	stopperigee.m_Condition = 0;  ///< 0=Cross Increasing;1=Cross Decreasing;2=Cross Either
	stopperigee.m_RepeatCount = 0;///< requested count, 期望重复次数
	stopperigee.m_TripValue = 0;  ///< 触发值(m, rad)
	stopperigee.m_Tolerance = 1e-6;///< 收敛控制误差(m, rad)
	double duration = 5184000;

	//调用高精度轨道外推，变步长，积分终止条件
	isSuccess = CalHPOPEphemeris(timePRL, false, PRL_LCOModOrbElem, 0, stopperigee, 2, duration, PosEJ2000, VelEJ2000, timeArray, posArray, velArray);
	time = timePRL;
	time.AddSec(duration);

	return isSuccess;
}


//***********************************************************
/// 由近月点伪状态正向求解返回出发时刻状态（双二体模型）
/// @Author	He Boyong   Zhou Wanmeng
/// @Date	2015-01-21  2018.5.12
/// @Input
/// @Param	timePRL						近月点时刻，UTCG
/// @Param	PRLHeight					近月点高度,m
/// @Param	PRLPerilunePseudoLat		近月点在月心瞬时惯性系中方位角纬度，rad
/// @Param	PRLPerilunePseudoLong		近月点在月心瞬时惯性系中方位角经度，rad
/// @Param	PRLPeriluneEccentric  		近月点LCO轨道偏心率，1.3，等价与近月点速度（m/s）一般在2.55e3 m/s
/// @Param	PRLPeriluneVelInc			近月点伪倾角，可以在pi±0.20943951，rad = 180±12，deg
/// @Output
/// @Param	LSOModElemEJ2000            飞出月球影响球时刻EJ2000的修正轨道根数
/// @Return		
//***********************************************************
void CFreeROInit::CalInitPseudoPRL2LSO(const CCTime& timePRL, const double& PRLHeight, const double& PRLPerilunePseudoLat, const double& PRLPerilunePseudoLong,
	const double& PRLPeriluneEccentric, double PRLPeriluneVelInc, CModOrbElem& LSOModElemEJ2000)
{
	CCTime tempTime, TT;							//近月点时刻
	CCoord3 posMJ2000, velMJ2000;				//近月点位置、速度
	CCoord3 posEJ2000, velEJ2000;				//出口点位置、速度
	CCoord3 lunarPosEJ2000, lunarVelEJ2000;     //月球的位置、速度
	double  timeLOI2LSO;						//飞出影响球的时间
	double  jdp;								//LSO时刻对应儒略日jd,

	CModOrbElem PRL_LCOModOrbElem;

	//由伪状态转换为近月点MJ2000位置速度
	PerilunePseudoState2MJ2000(timePRL, PRLHeight, PRLPerilunePseudoLat, PRLPerilunePseudoLong, PRLPeriluneEccentric, PRLPeriluneVelInc, posMJ2000, velMJ2000);
	//由位置速度计算近月点LCO修正轨道根数
	AsCartToModOrbElem(posMJ2000, velMJ2000, AsCMoonGrav, PRL_LCOModOrbElem);

	//参考时刻赋值
	tempTime = timePRL;

	//求解近月点到月球影响球出口点月心二体双曲线飞行时间
	CalTimeLSI2LOI(PRL_LCOModOrbElem, timeLOI2LSO);

	//获得LSO处的位置和速度矢量
	AsTwoBodyOrbProp(timeLOI2LSO, AsCMoonGrav, posMJ2000, velMJ2000);

	//从月心J2000坐标系月球影响球边界LSO转化为地心J2000坐标系参数
	tempTime.AddSec(timeLOI2LSO);

	//LSO时刻对应儒略日	
	AsUTCToTT(tempTime, TT);
	AsTimeToJD(TT, jdp);

	//JPL星历求解入口点瞬时月球相对于地球的位置速度
	AsJplPlanetaryEph(jdp, AsEJplMoon, AsEJplEarth, lunarPosEJ2000, lunarVelEJ2000);

	//月球影响球出口点地心J2000坐标系中状态
	posEJ2000 = posMJ2000 + lunarPosEJ2000;
	velEJ2000 = velMJ2000 + lunarVelEJ2000;

	//求解LSO处的EJ2000修正轨道根数
	AsCartToModOrbElem(posEJ2000, velEJ2000, AsCEarthGrav, LSOModElemEJ2000);

}
void CFreeROInit::CalInitElementsPRL2LSO(const CCTime& timePRL, const double& PRLHeight, const double& PRLPeriluneEccentric, const double& PRLPeriluneInc,
	const double& PRLPeriluneRAAN, double PRLPeriluneLAOP, CModOrbElem& LSOModElemEJ2000)
{
	CCTime tempTime, TT;							//近月点时刻
	CCoord3 posMJ2000, velMJ2000;				//近月点位置、速度
	CCoord3 posEJ2000, velEJ2000;				//出口点位置、速度
	CCoord3 lunarPosEJ2000, lunarVelEJ2000;     //月球的位置、速度
	double  timeLOI2LSO;						//飞出影响球的时间
	double  jdp;								//LSO时刻对应儒略日jd,

	CModOrbElem PRL_LCOModOrbElem;
	PRL_LCOModOrbElem.m_PeriRad = PRLHeight + AsCMoonRadius;
	PRL_LCOModOrbElem.m_Ecc = PRLPeriluneEccentric;
	PRL_LCOModOrbElem.m_I = PRLPeriluneInc;
	PRL_LCOModOrbElem.m_RAAN = PRLPeriluneRAAN;
	PRL_LCOModOrbElem.m_ArgPeri = PRLPeriluneLAOP;
	PRL_LCOModOrbElem.m_TrueA = 0;
	//参考时刻赋值
	tempTime = timePRL;

	//求解近月点到月球影响球出口点月心二体双曲线飞行时间
	CalTimeLSI2LOI(PRL_LCOModOrbElem, timeLOI2LSO);
	//由近月点LCO修正轨道根数计算位置速度
	AsModOrbElemToCart(PRL_LCOModOrbElem, AsCMoonGrav, posMJ2000, velMJ2000);
	//获得LSO处的位置和速度矢量
	AsTwoBodyOrbProp(timeLOI2LSO, AsCMoonGrav, posMJ2000, velMJ2000);

	//从月心J2000坐标系月球影响球边界LSO转化为地心J2000坐标系参数
	tempTime.AddSec(timeLOI2LSO);

	//LSO时刻对应儒略日	
	AsUTCToTT(tempTime, TT);
	AsTimeToJD(TT, jdp);

	//JPL星历求解入口点瞬时月球相对于地球的位置速度
	AsJplPlanetaryEph(jdp, AsEJplMoon, AsEJplEarth, lunarPosEJ2000, lunarVelEJ2000);

	//月球影响球出口点地心J2000坐标系中状态
	posEJ2000 = posMJ2000 + lunarPosEJ2000;
	velEJ2000 = velMJ2000 + lunarVelEJ2000;

	//求解LSO处的EJ2000修正轨道根数
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
	outfile.open(setFilepath.c_str());//与具体路径有关 需要改
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
		//近月点参数在地心J2000坐标系的参数
		AsJplPlanetaryEph(jd, AsEJplMoon, AsEJplEarth, lunarPosEJ2000, lunarVelEJ2000);
		//近月点相对地球的EJ2000位置速度
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

		///高精度结果
		outfile
			<< "============================" << endl
			<< "自由返回轨道的仿真结果报告" << i << endl
			<< "============================" << endl << endl
			<< "发射转移时间（sec）:   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.backwardTime[i] << "   " << endl << endl
			<< "发射转移地心夹角（deg）:   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLIPhaseError[i] * AsCRadToDeg << "   " << endl << endl
			<< "-----------" << endl
			<< "入射点状态" << endl
			<< "-----------" << endl << endl
			<< "初始历元(UTCG):   " << endl << endl
			<< setw(18) << setprecision(15) << launchTimestr << "   " << endl << endl

			<< "初始历元(UTCG):   " << endl << endl
			<< setw(18) << setprecision(15) << launchTimestr << "   " << endl << endl
			<< "变轨脉冲 VVLH（m/s）:   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLIDeltaV[i] << "   "
			<< setw(18) << setprecision(15) << 0 << "  "
			<< setw(18) << setprecision(15) << 0 << endl << endl
			<< "EJ2000位置(m)、速度(m/s):   " << endl << endl
			<< "x:" << "  " << setw(18) << setprecision(15) << PosEJ2000[0] << endl
			<< "y:" << "  " << setw(18) << setprecision(15) << PosEJ2000[1] << endl
			<< "z:" << "  " << setw(18) << setprecision(15) << PosEJ2000[2] << endl
			<< "vx:" << "  " << setw(18) << setprecision(15) << VelEJ2000[0] << endl
			<< "vy:" << "  " << setw(18) << setprecision(15) << VelEJ2000[1] << endl
			<< "vz:" << "  " << setw(18) << setprecision(15) << VelEJ2000[2] << endl << endl
			<< "EJ2000修正轨道根数:   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_PeriRad << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_Ecc << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_I*AsCRadToDeg << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_RAAN*AsCRadToDeg << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_ArgPeri*AsCRadToDeg << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_TrueA*AsCRadToDeg << endl << endl
			<< "EJ2000轨道根数:   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_PeriRad / (1 - mFreeRTONRCalOutput.TLI_LTOModElem[i].m_Ecc) << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_Ecc << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_I*AsCRadToDeg << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_RAAN*AsCRadToDeg << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_ArgPeri*AsCRadToDeg << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TLI_LTOModElem[i].m_TrueA*AsCRadToDeg << endl << endl
			<< "----------" << endl
			<< "近月点状态" << endl
			<< "----------" << endl << endl
			<< "初始历元(UTCG):   " << endl << endl
			<< setw(18) << setprecision(15) << PRLtimestr << endl << endl
			<< "EJ2000位置(m)、速度(m/s):   " << endl << endl
			<< "x:" << setw(18) << setprecision(15) << PRL_LCOPosEJ2000[0] << endl
			<< "y:" << setw(18) << setprecision(15) << PRL_LCOPosEJ2000[1] << endl
			<< "z:" << setw(18) << setprecision(15) << PRL_LCOPosEJ2000[2] << endl
			<< "vx:" << setw(18) << setprecision(15) << PRL_LCOVelEJ2000[0] << endl
			<< "vy:" << setw(18) << setprecision(15) << PRL_LCOVelEJ2000[1] << endl
			<< "vz:" << setw(18) << setprecision(15) << PRL_LCOVelEJ2000[2] << endl << endl
			<< "MJ2000修正轨道根数:   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LCOModElem[i].m_PeriRad << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LCOModElem[i].m_Ecc << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LCOModElem[i].m_I *AsCRadToDeg << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LCOModElem[i].m_RAAN *AsCRadToDeg << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LCOModElem[i].m_ArgPeri *AsCRadToDeg << "   "
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LCOModElem[i].m_TrueA *AsCRadToDeg << endl << endl
			<< "----------" << endl
			<< "再入点状态" << endl
			<< "----------" << endl << endl
			<< "返回轨道倾角(deg):   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.RTOInc[i] * AsCRadToDeg << endl << endl
			<< "返回点经度差(deg):   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.ReturnLonError[i] * AsCRadToDeg << endl << endl
			<< "真空近地点高度（m）:   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.VCPH[i] << endl << endl
			<< "再入历元(UTCG):   " << endl << endl
			<< setw(18) << setprecision(15) << ROItimestr << endl << endl
			<< "再入角(deg):   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.ReentryGamma[i] * AsCRadToDeg << endl << endl
			<< "再入点地理纬度(deg):   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.ROIECFLat[i] * AsCRadToDeg << endl << endl
			<< "再入点地理经度(deg):   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.ROIECFLon[i] * AsCRadToDeg << endl << endl
			<< "再入点地心角(deg):   " << endl << endl
			<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.ROIPhaseError[i] * AsCRadToDeg << endl << endl;
		if (!mIsImpulse)
		{
			outfile << endl;
		}
		else
		{
			outfile
				<< "----------" << endl
				<< "第一次变轨" << endl
				<< "----------" << endl << endl
				<< "变轨历元（UTCG）:   " << endl << endl
				<< setw(18) << setprecision(15) << timestr1 << endl << endl
				<< "变轨脉冲 MJ2000（m/s）:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOIDeltaV[i][0] << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOIDeltaV[i][1] << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOIDeltaV[i][2] << endl << endl
				<< "转移前MJ2000修正轨道根数:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOBeforeModElem[i].m_PeriRad << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOBeforeModElem[i].m_Ecc << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOBeforeModElem[i].m_I* AsCRadToDeg << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOBeforeModElem[i].m_RAAN * AsCRadToDeg << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOBeforeModElem[i].m_ArgPeri* AsCRadToDeg << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOBeforeModElem[i].m_TrueA* AsCRadToDeg << endl << endl
				<< "MJ2000修正轨道根数:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOModElem[i].m_PeriRad << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOModElem[i].m_Ecc << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOModElem[i].m_I *AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOModElem[i].m_RAAN *AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOModElem[i].m_ArgPeri *AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.PRL_LIOModElem[i].m_TrueA *AsCRadToDeg << endl << endl
				<< "----------" << endl
				<< "第二次变轨" << endl
				<< "----------" << endl << endl
				<< "变轨历元（UTCG）:   " << endl << endl
				<< setw(18) << setprecision(15) << timestr2 << endl << endl
				<< "变轨脉冲 MJ2000（m/s）:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2DeltaV[i][0] << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2DeltaV[i][1] << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2DeltaV[i][2] << endl << endl
				<< "转移前MJ2000修正轨道根数:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2BeforeModElem[i].m_PeriRad << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2BeforeModElem[i].m_Ecc << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2BeforeModElem[i].m_I* AsCRadToDeg << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2BeforeModElem[i].m_RAAN * AsCRadToDeg << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2BeforeModElem[i].m_ArgPeri* AsCRadToDeg << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2BeforeModElem[i].m_TrueA* AsCRadToDeg << endl << endl
				<< "MJ2000修正轨道根数:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2ModElem[i].m_PeriRad << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2ModElem[i].m_Ecc << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2ModElem[i].m_I *AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2ModElem[i].m_RAAN *AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2ModElem[i].m_ArgPeri *AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI2ModElem[i].m_TrueA *AsCRadToDeg << endl << endl
				<< "----------" << endl
				<< "第三次变轨" << endl
				<< "----------" << endl << endl
				<< "变轨历元（UTCG）:   " << endl << endl
				<< setw(18) << setprecision(15) << timestr3 << endl << endl
				<< "变轨脉冲 MJ2000（m/s）:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3DeltaV[i][0] << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3DeltaV[i][1] << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3DeltaV[i][2] << endl << endl
				<< "转移前MJ2000修正轨道根数:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3BeforeModElem[i].m_PeriRad << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3BeforeModElem[i].m_Ecc << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3BeforeModElem[i].m_I* AsCRadToDeg << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3BeforeModElem[i].m_RAAN * AsCRadToDeg << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3BeforeModElem[i].m_ArgPeri* AsCRadToDeg << "  "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3BeforeModElem[i].m_TrueA* AsCRadToDeg << endl << endl
				<< "MJ2000修正轨道根数:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3ModElem[i].m_PeriRad << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3ModElem[i].m_Ecc << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3ModElem[i].m_I *AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3ModElem[i].m_RAAN *AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3ModElem[i].m_ArgPeri *AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LOI3ModElem[i].m_TrueA *AsCRadToDeg << endl << endl
				<< "------------" << endl
				<< "变轨后的状态" << endl
				<< "------------" << endl << endl
				<< "结束历元（UTCG）:   " << endl << endl
				<< setw(18) << setprecision(15) << timeend << endl << endl
				<< "MJ2000修正轨道根数:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LDOModElem[i].m_PeriRad << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LDOModElem[i].m_Ecc << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LDOModElem[i].m_I*AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LDOModElem[i].m_RAAN*AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LDOModElem[i].m_ArgPeri*AsCRadToDeg << "   "
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.LDOModElem[i].m_TrueA*AsCRadToDeg << endl << endl
				<< "总速度增量（m/s）:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.TotalImpulse[i] << endl << endl
				<< "异面差（deg）:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.epsi[i] * AsCRadToDeg << endl << endl
				<< "拱线夹角（deg）:   " << endl << endl
				<< setw(18) << setprecision(15) << mFreeRTONRCalOutput.perierror[i] * AsCRadToDeg << endl << endl;
		}
	}
	outfile.close();
	return true;
}


















