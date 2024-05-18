#pragma once
#include"MannedLunarCommon.h"
#include"Impulses.h"
#include "AstroLib.h"
#include<string>



////////////////////////////////////////////////////////////////////
///////////////////自由返回转移轨道（FreeReturnOrbit）///////////////////////
////////////////////////////////////////////////////////////////////

//初步计算输出参数
struct FreeReturnOrbitInitCalOutput
{
	std::vector<CCTime> TLILaunchTime;      //TLI时刻（UTC）
	std::vector<double> TLIDeltaV;		//<<TLI入射点脉冲速度绝对值（m/s）
	std::vector<CModOrbElem>TLI_LTOModElem;//近地出发轨道根数
	std::vector<double> LTOInc;             //LTO倾角（rad）
	std::vector<CCTime> PRLTime;            //近月点时刻
	std::vector<CModOrbElem>LCO_PRLModElem; //近月制动前的轨道根数
	std::vector<CModOrbElem>LSIModElem;     //入口点前的轨道根数
	std::vector<CModOrbElem>LSOModElem;     //
	std::vector<double>PRLPseudoLat;        //近月点时刻近月点伪纬度
	std::vector<double>PRLPseudoLong;       //近月点时刻近月段伪经度
	std::vector<double>PRLVelInc;           //近月点时刻近月点速度倾角
	std::vector<double>PRLEccentric;        //近月点时刻近月点轨道偏心率
	std::vector<bool>IsLTODesend;           //LTO轨道升轨还是降轨
	std::vector<bool>IsLDODesend;           //月心段轨道升轨还是降轨
	std::vector<bool>IsRTODesend;           //返回轨道升轨还是降轨
	std::vector<double>TimeofFlight;        //飞行总时间
	std::vector<double>ROIInc;              //返回轨道倾角
	std::vector<double>VCPH;                //真空近地点高度
	std::vector<CCTime>ROITime;		    	//ROI时刻 (UTC)
	std::vector<double>ROILat;              //ROI纬度
	std::vector<double>ROILon;              //ROI经度

};
//高精度计算输出参数
struct NRHPFreeRTOCalOutput
{
	std::vector<CCTime>      TLILaunchTime;           ////<<TLI时刻 (UTC)
	std::vector<double>      TLIDeltaV;				 ////<<TLI入射点脉冲速度绝对值（m/s）
	std::vector<CModOrbElem> TLI_LTOModElem;          ////<<TLI时刻，LTO的轨道修正根数，详见AstroLib说明
	std::vector<CModOrbElem> TLI_LTOModElemE;         ////<<TLI时刻，地固系LTO的轨道修正根数，详见AstroLib说明

	std::vector<CCTime>      PRLTime;		    	 ////<<PRL/LOI时刻 (UTC)
	std::vector<CModOrbElem> PRL_LCOModElem;         ////<<PRL时刻，LCO的轨道修正根数，详见AstroLib说明
	std::vector<CModOrbElem> LLOModElem;             ////<<PRL时刻，LLO的轨道修正根数，详见AstroLib说明
	std::vector<double>      RTOInc;				 ////<<RTO倾角(rad) 
	std::vector<double>      VCPH;                   ////<<真空近地点高度（m）
	std::vector<double>      epsi;
	std::vector<double>      perierror;
	std::vector<double>      ReentryGamma;

	std::vector<CCoord3>	 LOIDeltaV;				 ////<<环月1脉冲捕获制动速度 j2000（m/s）
	std::vector<CModOrbElem> PRL_LIOModElem;
	std::vector<CModOrbElem> PRL_LIOBeforeModElem;
	std::vector<CCTime>		 LOI1Time;		    	 ////<<LOI1时刻 (UTC)
	std::vector<CCTime>		 LOI2Time;		    	 ////<<LOI2时刻 (UTC)
	std::vector<CCoord3>	 LOI2DeltaV;			 ////<<环月2脉冲捕获制动速度 j2000（m/s）
	std::vector<CModOrbElem> LOI2ModElem;            ////<<二脉冲后，LIO的轨道修正根数，详见AstroLib说明
	std::vector<CModOrbElem> LOI2BeforeModElem;
	std::vector<CCTime>		 LOI3Time;		    	 ////<<LOI3时刻 (UTC)
	std::vector<CCoord3>	 LOI3DeltaV;		     ////<<环月3脉冲捕获制动速度 j2000（m/s）
	std::vector<CModOrbElem> LOI3ModElem;
	std::vector<CModOrbElem> LOI3BeforeModElem;
	std::vector<CModOrbElem> LDOModElem;             ////<<RVD开始时刻，LOI3的轨道修正根数，详见AstroLib说明

	std::vector<double>		 TotalImpulse;           ////<<3脉冲转移的总速度增量（m/s）

	std::vector<CCTime>		 ROITime;		    	////<<ROI时刻 (UTC)
	std::vector<double>		 ROIECFLat;             ////<<ROI,点地理纬度（rad）
	std::vector<double>		 ROIECFLon;             ////<<ROI,点地理经度（rad）
	std::vector<double>		 ROIPhaseError;         ////<<ROI距着陆点地心夹角（rad）
	std::vector<double>      ReturnLonError;        ////<<返回场经度的差（rad）

	std::vector<double>      TLIPhaseError;
	std::vector<double>      backwardTime;
};
///******自由返回转移轨道（FreeReturnOrbit）******///
class CFreeROInit
{
public:
	CFreeROInit();
	~CFreeROInit();

	//初步计算成员变量
	//时间设计参数
	CCTime mPRLTime[2];           //近月点时刻
	double mStep;                 //遍历步长
	double mTof;                  //地月空间转移时间
	double mPRL2RVDTime;
	bool mIsLCODesend;            //月心段升降轨
	double mLCO_PRLHeight;        //近月点高度

	//TLI时刻LTO轨道参数
	double mLTOPeriHeight;        //LTO近地点高度
	double mLTOInc;               //LTO轨道倾角
	bool mIsLTODesend;            //升降轨发射
	double mLTOLAN;               //升交点经度
	double mLTOTrueA;             //真近点角

	//ROI时刻，RTO轨道参数
	double mReentryGamma;         //再入角
	double mAtmospH;              //大气层高度
	double mRTOInc;               //返回轨道倾角
	bool mIsRTODesend;            //RTO轨道升降轨

	double mLaunchSiteLat;        //发射场纬度
	double mLaunchSiteLon;        //发射场经度
	double mLaunchHeight;         //发射场高度
	double mReturnSiteECFLon;     //着陆场地理经度
	double mReturnSiteECFLat;     //着陆场地理纬度

	///计算输出参数
	FreeReturnOrbitInitCalOutput mFreeReturnOrbitInitCalOutput;  //初步轨道计算结果
	//PRL时刻LDO(MJ2000)参数
	CCTime       mLOITime;                 ///<<近月点时刻
	CCTime       mTLITime;                 ///<<出发时刻

	//来自年窗口计算的着陆时刻和环月LLO轨道参数
	CCTime       mLLOTime;                 ///<<参考环月轨道历元（UTC）
	CModOrbElem  mLLOModElem;              ///<<环月LLO轨道参数，详见AstroLib说明
	CModOrbElem  mLLOModElemNew;           ///<<环月LLO轨道参数，详见AstroLib说明
	//PRL时刻，设计参数
	double mPRLPseudoLat;         //近月点伪纬度
	double mPRLPseudoLon;         //近月点伪经度
	double mPRLVelInc;            //近月点速度倾角
	double mPRLEccentric;         //近月点轨道偏心率

	double PRLPeriluneEccentric;         //月心段轨道偏心率
	double PRLPeriluneInc;         //月心段轨道倾角
	double PRLPeriluneRAAN;            //月心段轨道升交点赤经
	double PRLPeriluneLAOP;         //月心段轨道近拱点纬度幅角

	double mLTOPeriHeightError;
	double mLTOIncError;
	double mRTOIncError;
	double mVCPHError;
	double mTimeError;

	double mVCPH;
	double mTrueAnalROI;
	double mTofTemp;              //地月空间转移时间

	CModOrbElem mLSIModElemEJ2000;
	CModOrbElem mLSOModElemEJ2000;
	CModOrbElem mLCO_PRLModElem;
	bool mIsSuccess;
	bool mIsImpulse;
	double mHErrorLimit;           //LTO近地点高度误差限，10m
	double mIncErrorLimit;         //LTO轨道倾角
	double mTOFLimit;              //TOF不超过3天

	CCTime mROITime;               //再入点时刻
	double mROIECFLat;             //再入点地理纬度
	double mROIECFLon;             //再入点地理经度

	double mepsi, mperierror;
	int mloop;
	//
	//RVD时刻LDO(MJ2000)参数
	double	     mLDOEcc;		           ///<偏心率 eccentricity
	double	     mLDOArgPeri;	           ///<近拱点角 arg of periapsis in rad
	double	     mLDOTrueA;	               ///<真近点角 true anomaly in radians
	double       malpha;                   ///<中间转移轨道的周期比例

	//ROI时刻，RTO轨道，约束参数
	double       mVCPDecli;                ///<<真空近地点赤纬（rad）
	double       mThrustTime;                ///<<发动机工作时间（sec）,1500 s
	//深空站位置
	std::vector<CLatLonAlt> mFacLLA;
	std::vector<double>     mAngConstr;

	///计算输出参数
	NRHPFreeRTOCalOutput   mFreeRTONRCalOutput;     ///<<高精度轨道计算结果
	std::string           mResultFilepath;	        ///输出参数文件路径
	///格式: 日期
	///格式：TLI时刻 (UTC)，TLI deltaV,TLI时刻LTOModElem,LOI时刻 (UTC),LOI deltaV，
	///	     PRL时刻近月点纬度,PRL时刻近月点经度,PRL时刻近月点速度倾角,PRL时刻近月点轨道偏心率
	///      LOI时刻后LIOModElem,LOI2时刻 (UTC),LOI2 deltaV,LOI3时刻 (UTC),LOI3 deltaV，
	///		 RVD开始时刻LIOModElem，轨道类型
	///曲线数据
	std::vector<double>     mTime;              ///<<遍历时间（s）
	std::vector<CCoord3> 	mEJ2000Pos;   		///<<地心J2000（m）
	std::vector<CCoord3> 	mEJ2000Vel;   		///<<地心J2000（m/s）
	std::vector<CCoord3> 	mMJ2000Pos;   		///<<月心J2000（m）
	std::vector<CCoord3> 	mMJ2000Vel;   		///<<月心J2000（m/s）
	std::vector<bool>       mIsCovered;

	CCoord3			mLOIDeltaV;				 ////<<环月1脉冲捕获制动速度 j2000（m/s）
	CModOrbElem		mPRL_LIOModElem;
	CCTime			mLOI2Time;		    	 ////<<LOI2时刻 (UTC)
	CCoord3			mLOI2DeltaV;			 ////<<环月2脉冲捕获制动速度 j2000（m/s）
	CModOrbElem		mLOI2ModElem;            ////<<二脉冲后，LIO的轨道修正根数，详见AstroLib说明

	CCTime			mLOI3Time;		    	 ////<<LOI3时刻 (UTC)
	CCoord3			mLOI3DeltaV;		     ////<<环月3脉冲捕获制动速度 j2000（m/s）
	CModOrbElem		mLOI3ModElem;            ////<<RVD开始时刻，LOI3的轨道修正根数，详见AstroLib说明
	CModOrbElem		mLDOModElem;
	double			mTotalImpulse;           ////<<3脉冲转移的总速度增量（m/s）



	CModOrbElem mTLIModOrbElem;
	CModOrbElem mVCPModOrbElem;


	CVector<double> int_interval;
	CVector<double> ks;
	CVector<double> epsino;
	CVector<double> Cg_Optx;
	//成员函数
	typedef void(CFreeROInit::*Senfunc)(double&, double&);
	typedef void(CFreeROInit::*Calfunc)(CVector<double>&, CVector<double>&);
	typedef void(CFreeROInit::*Confunc)(CVector<double>, CVector<double>, CVector<double>, bool&);
	typedef void(CFreeROInit::*SenMa)(Calfunc,  CVector<double>&, CVector<double>&,
		CVector<double>, CVector<double>,  CMatrix<double>&);
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
	void CalInitPseudoPRL2LSI(const CCTime& timePRL, const double& PRLHeight, const double& PRLPerilunePseudoLat, const double& PRLPerilunePseudoLong,
		const double& PRLPeriluneEccentric, double PRLPeriluneVelInc, CModOrbElem& LSIModElemEJ2000);
	void CalInitElementsPRL2LSI(const CCTime& timePRL, const double& PRLHeight, const double& PRLPeriluneEccentric, const double& PRLPeriluneInc,
		const double& PRLPeriluneRAAN, double PRLPeriluneLAOP, CModOrbElem& LSIModElemEJ2000);
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
	void PerilunePseudoState2MJ2000(const CCTime& timePRL, const double& PRLHeight, const double& PRLPerilunePseudoLat,
		const double& PRLPerilunePseudoLong, const double& PRLPeriluneEccentric, const double& PRLPeriluneVelInc,
		CCoord3& posMJ2000, CCoord3& velMJ2000);
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
	/// @Param	LSOModElemEJ2000      进入月球影响球时刻EJ2000的修正轨道根数
	/// @Return		
	//***********************************************************
	void CalInitPseudoPRL2LSO(const CCTime& timePRL, const double& PRLHeight, const double& PRLPerilunePseudoLat, const double& PRLPerilunePseudoLong,
		const double& PRLPeriluneEccentric, double PRLPeriluneVelInc, CModOrbElem& LSOModElemEJ2000);
	void CalInitElementsPRL2LSO(const CCTime& timePRL, const double& PRLHeight, const double& PRLPeriluneEccentric, const double& PRLPeriluneInc,
		const double& PRLPeriluneRAAN, double PRLPeriluneLAOP, CModOrbElem& LSOModElemEJ2000);
	//***********************************************************
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
	bool CalHPLOI2VCP(const CCTime& timePRL, const CCoord3& PRL_LCOPosMJ2000, const CCoord3& PRL_LCOVelMJ2000,
		CCTime& time, CCoord3& PosEJ2000, CCoord3& VelEJ2000, std::vector<double>& timeArray,
		std::vector<CCoord3>& posArray, std::vector<CCoord3>& velArray);

	void SenCalFun( CVector<double>& x, CVector<double>& result);
	void ConJud(CVector<double> Out_Y, CVector<double> Y_tar, CVector<double> Out_Ylim, bool&IsCon);
	void SenCal(Calfunc f,  CVector<double>& x, CVector<double>& result, 
		CVector<double> int_interval, CVector<double> k, CMatrix<double>& sen_result);
	void NRMethod(Confunc f1, Calfunc f2, SenMa f3, Calfunc f4, CVector<double>& x,
		CVector<double> Y_tar, CVector<double> Out_Ylim, CVector<double> int_interval,
		CMatrix<double> Iter_interval, CVector<double>k,int Iternum, bool&IsCon);

	void NRFreeRTOInit(const CModOrbElem& LLOModOrbElem, CVector<double>& x,
		CVector<double> Y_tar, CVector<double> Out_Ylim, CVector<double> int_interval,
		CMatrix<double> Iter_interval, CVector<double> k, int Iternum, bool&IsCon);
	void NRFreeRTOHP(const CModOrbElem& LLOModOrbElem, CVector<double>& x,
		CVector<double> Y_tar, CVector<double> Out_Ylim, CVector<double> int_interval,
		CMatrix<double> Iter_interval, CVector<double> k, int Iternum, bool&IsCon);
	void SenHPCalFun(CVector<double>& x, CVector<double>& result);
	bool HPSaveFreeRTOResult(const std::string& setFilepath);
	void SenCalFunElements(CVector<double>& x, CVector<double>& result);
	void SenHPCalFunElement(CVector<double>& x, CVector<double>& result);
};


