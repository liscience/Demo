#pragma once
#include"MannedLunarCommon.h"
#include"Impulses.h"
#include "AstroLib.h"
#include<string>



////////////////////////////////////////////////////////////////////
///////////////////���ɷ���ת�ƹ����FreeReturnOrbit��///////////////////////
////////////////////////////////////////////////////////////////////

//���������������
struct FreeReturnOrbitInitCalOutput
{
	std::vector<CCTime> TLILaunchTime;      //TLIʱ�̣�UTC��
	std::vector<double> TLIDeltaV;		//<<TLI����������ٶȾ���ֵ��m/s��
	std::vector<CModOrbElem>TLI_LTOModElem;//���س����������
	std::vector<double> LTOInc;             //LTO��ǣ�rad��
	std::vector<CCTime> PRLTime;            //���µ�ʱ��
	std::vector<CModOrbElem>LCO_PRLModElem; //�����ƶ�ǰ�Ĺ������
	std::vector<CModOrbElem>LSIModElem;     //��ڵ�ǰ�Ĺ������
	std::vector<CModOrbElem>LSOModElem;     //
	std::vector<double>PRLPseudoLat;        //���µ�ʱ�̽��µ�αγ��
	std::vector<double>PRLPseudoLong;       //���µ�ʱ�̽��¶�α����
	std::vector<double>PRLVelInc;           //���µ�ʱ�̽��µ��ٶ����
	std::vector<double>PRLEccentric;        //���µ�ʱ�̽��µ���ƫ����
	std::vector<bool>IsLTODesend;           //LTO������컹�ǽ���
	std::vector<bool>IsLDODesend;           //���Ķι�����컹�ǽ���
	std::vector<bool>IsRTODesend;           //���ع�����컹�ǽ���
	std::vector<double>TimeofFlight;        //������ʱ��
	std::vector<double>ROIInc;              //���ع�����
	std::vector<double>VCPH;                //��ս��ص�߶�
	std::vector<CCTime>ROITime;		    	//ROIʱ�� (UTC)
	std::vector<double>ROILat;              //ROIγ��
	std::vector<double>ROILon;              //ROI����

};
//�߾��ȼ����������
struct NRHPFreeRTOCalOutput
{
	std::vector<CCTime>      TLILaunchTime;           ////<<TLIʱ�� (UTC)
	std::vector<double>      TLIDeltaV;				 ////<<TLI����������ٶȾ���ֵ��m/s��
	std::vector<CModOrbElem> TLI_LTOModElem;          ////<<TLIʱ�̣�LTO�Ĺ���������������AstroLib˵��
	std::vector<CModOrbElem> TLI_LTOModElemE;         ////<<TLIʱ�̣��ع�ϵLTO�Ĺ���������������AstroLib˵��

	std::vector<CCTime>      PRLTime;		    	 ////<<PRL/LOIʱ�� (UTC)
	std::vector<CModOrbElem> PRL_LCOModElem;         ////<<PRLʱ�̣�LCO�Ĺ���������������AstroLib˵��
	std::vector<CModOrbElem> LLOModElem;             ////<<PRLʱ�̣�LLO�Ĺ���������������AstroLib˵��
	std::vector<double>      RTOInc;				 ////<<RTO���(rad) 
	std::vector<double>      VCPH;                   ////<<��ս��ص�߶ȣ�m��
	std::vector<double>      epsi;
	std::vector<double>      perierror;
	std::vector<double>      ReentryGamma;

	std::vector<CCoord3>	 LOIDeltaV;				 ////<<����1���岶���ƶ��ٶ� j2000��m/s��
	std::vector<CModOrbElem> PRL_LIOModElem;
	std::vector<CModOrbElem> PRL_LIOBeforeModElem;
	std::vector<CCTime>		 LOI1Time;		    	 ////<<LOI1ʱ�� (UTC)
	std::vector<CCTime>		 LOI2Time;		    	 ////<<LOI2ʱ�� (UTC)
	std::vector<CCoord3>	 LOI2DeltaV;			 ////<<����2���岶���ƶ��ٶ� j2000��m/s��
	std::vector<CModOrbElem> LOI2ModElem;            ////<<�������LIO�Ĺ���������������AstroLib˵��
	std::vector<CModOrbElem> LOI2BeforeModElem;
	std::vector<CCTime>		 LOI3Time;		    	 ////<<LOI3ʱ�� (UTC)
	std::vector<CCoord3>	 LOI3DeltaV;		     ////<<����3���岶���ƶ��ٶ� j2000��m/s��
	std::vector<CModOrbElem> LOI3ModElem;
	std::vector<CModOrbElem> LOI3BeforeModElem;
	std::vector<CModOrbElem> LDOModElem;             ////<<RVD��ʼʱ�̣�LOI3�Ĺ���������������AstroLib˵��

	std::vector<double>		 TotalImpulse;           ////<<3����ת�Ƶ����ٶ�������m/s��

	std::vector<CCTime>		 ROITime;		    	////<<ROIʱ�� (UTC)
	std::vector<double>		 ROIECFLat;             ////<<ROI,�����γ�ȣ�rad��
	std::vector<double>		 ROIECFLon;             ////<<ROI,������ȣ�rad��
	std::vector<double>		 ROIPhaseError;         ////<<ROI����½����ļнǣ�rad��
	std::vector<double>      ReturnLonError;        ////<<���س����ȵĲrad��

	std::vector<double>      TLIPhaseError;
	std::vector<double>      backwardTime;
};
///******���ɷ���ת�ƹ����FreeReturnOrbit��******///
class CFreeROInit
{
public:
	CFreeROInit();
	~CFreeROInit();

	//���������Ա����
	//ʱ����Ʋ���
	CCTime mPRLTime[2];           //���µ�ʱ��
	double mStep;                 //��������
	double mTof;                  //���¿ռ�ת��ʱ��
	double mPRL2RVDTime;
	bool mIsLCODesend;            //���Ķ�������
	double mLCO_PRLHeight;        //���µ�߶�

	//TLIʱ��LTO�������
	double mLTOPeriHeight;        //LTO���ص�߶�
	double mLTOInc;               //LTO������
	bool mIsLTODesend;            //�����췢��
	double mLTOLAN;               //�����㾭��
	double mLTOTrueA;             //������

	//ROIʱ�̣�RTO�������
	double mReentryGamma;         //�����
	double mAtmospH;              //������߶�
	double mRTOInc;               //���ع�����
	bool mIsRTODesend;            //RTO���������

	double mLaunchSiteLat;        //���䳡γ��
	double mLaunchSiteLon;        //���䳡����
	double mLaunchHeight;         //���䳡�߶�
	double mReturnSiteECFLon;     //��½��������
	double mReturnSiteECFLat;     //��½������γ��

	///�����������
	FreeReturnOrbitInitCalOutput mFreeReturnOrbitInitCalOutput;  //�������������
	//PRLʱ��LDO(MJ2000)����
	CCTime       mLOITime;                 ///<<���µ�ʱ��
	CCTime       mTLITime;                 ///<<����ʱ��

	//�����괰�ڼ������½ʱ�̺ͻ���LLO�������
	CCTime       mLLOTime;                 ///<<�ο����¹����Ԫ��UTC��
	CModOrbElem  mLLOModElem;              ///<<����LLO������������AstroLib˵��
	CModOrbElem  mLLOModElemNew;           ///<<����LLO������������AstroLib˵��
	//PRLʱ�̣���Ʋ���
	double mPRLPseudoLat;         //���µ�αγ��
	double mPRLPseudoLon;         //���µ�α����
	double mPRLVelInc;            //���µ��ٶ����
	double mPRLEccentric;         //���µ���ƫ����

	double PRLPeriluneEccentric;         //���Ķι��ƫ����
	double PRLPeriluneInc;         //���Ķι�����
	double PRLPeriluneRAAN;            //���Ķι��������ྭ
	double PRLPeriluneLAOP;         //���Ķι��������γ�ȷ���

	double mLTOPeriHeightError;
	double mLTOIncError;
	double mRTOIncError;
	double mVCPHError;
	double mTimeError;

	double mVCPH;
	double mTrueAnalROI;
	double mTofTemp;              //���¿ռ�ת��ʱ��

	CModOrbElem mLSIModElemEJ2000;
	CModOrbElem mLSOModElemEJ2000;
	CModOrbElem mLCO_PRLModElem;
	bool mIsSuccess;
	bool mIsImpulse;
	double mHErrorLimit;           //LTO���ص�߶�����ޣ�10m
	double mIncErrorLimit;         //LTO������
	double mTOFLimit;              //TOF������3��

	CCTime mROITime;               //�����ʱ��
	double mROIECFLat;             //��������γ��
	double mROIECFLon;             //����������

	double mepsi, mperierror;
	int mloop;
	//
	//RVDʱ��LDO(MJ2000)����
	double	     mLDOEcc;		           ///<ƫ���� eccentricity
	double	     mLDOArgPeri;	           ///<������� arg of periapsis in rad
	double	     mLDOTrueA;	               ///<������ true anomaly in radians
	double       malpha;                   ///<�м�ת�ƹ�������ڱ���

	//ROIʱ�̣�RTO�����Լ������
	double       mVCPDecli;                ///<<��ս��ص��γ��rad��
	double       mThrustTime;                ///<<����������ʱ�䣨sec��,1500 s
	//���վλ��
	std::vector<CLatLonAlt> mFacLLA;
	std::vector<double>     mAngConstr;

	///�����������
	NRHPFreeRTOCalOutput   mFreeRTONRCalOutput;     ///<<�߾��ȹ��������
	std::string           mResultFilepath;	        ///��������ļ�·��
	///��ʽ: ����
	///��ʽ��TLIʱ�� (UTC)��TLI deltaV,TLIʱ��LTOModElem,LOIʱ�� (UTC),LOI deltaV��
	///	     PRLʱ�̽��µ�γ��,PRLʱ�̽��µ㾭��,PRLʱ�̽��µ��ٶ����,PRLʱ�̽��µ���ƫ����
	///      LOIʱ�̺�LIOModElem,LOI2ʱ�� (UTC),LOI2 deltaV,LOI3ʱ�� (UTC),LOI3 deltaV��
	///		 RVD��ʼʱ��LIOModElem���������
	///��������
	std::vector<double>     mTime;              ///<<����ʱ�䣨s��
	std::vector<CCoord3> 	mEJ2000Pos;   		///<<����J2000��m��
	std::vector<CCoord3> 	mEJ2000Vel;   		///<<����J2000��m/s��
	std::vector<CCoord3> 	mMJ2000Pos;   		///<<����J2000��m��
	std::vector<CCoord3> 	mMJ2000Vel;   		///<<����J2000��m/s��
	std::vector<bool>       mIsCovered;

	CCoord3			mLOIDeltaV;				 ////<<����1���岶���ƶ��ٶ� j2000��m/s��
	CModOrbElem		mPRL_LIOModElem;
	CCTime			mLOI2Time;		    	 ////<<LOI2ʱ�� (UTC)
	CCoord3			mLOI2DeltaV;			 ////<<����2���岶���ƶ��ٶ� j2000��m/s��
	CModOrbElem		mLOI2ModElem;            ////<<�������LIO�Ĺ���������������AstroLib˵��

	CCTime			mLOI3Time;		    	 ////<<LOI3ʱ�� (UTC)
	CCoord3			mLOI3DeltaV;		     ////<<����3���岶���ƶ��ٶ� j2000��m/s��
	CModOrbElem		mLOI3ModElem;            ////<<RVD��ʼʱ�̣�LOI3�Ĺ���������������AstroLib˵��
	CModOrbElem		mLDOModElem;
	double			mTotalImpulse;           ////<<3����ת�Ƶ����ٶ�������m/s��



	CModOrbElem mTLIModOrbElem;
	CModOrbElem mVCPModOrbElem;


	CVector<double> int_interval;
	CVector<double> ks;
	CVector<double> epsino;
	CVector<double> Cg_Optx;
	//��Ա����
	typedef void(CFreeROInit::*Senfunc)(double&, double&);
	typedef void(CFreeROInit::*Calfunc)(CVector<double>&, CVector<double>&);
	typedef void(CFreeROInit::*Confunc)(CVector<double>, CVector<double>, CVector<double>, bool&);
	typedef void(CFreeROInit::*SenMa)(Calfunc,  CVector<double>&, CVector<double>&,
		CVector<double>, CVector<double>,  CMatrix<double>&);
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
	void CalInitPseudoPRL2LSI(const CCTime& timePRL, const double& PRLHeight, const double& PRLPerilunePseudoLat, const double& PRLPerilunePseudoLong,
		const double& PRLPeriluneEccentric, double PRLPeriluneVelInc, CModOrbElem& LSIModElemEJ2000);
	void CalInitElementsPRL2LSI(const CCTime& timePRL, const double& PRLHeight, const double& PRLPeriluneEccentric, const double& PRLPeriluneInc,
		const double& PRLPeriluneRAAN, double PRLPeriluneLAOP, CModOrbElem& LSIModElemEJ2000);
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
	void PerilunePseudoState2MJ2000(const CCTime& timePRL, const double& PRLHeight, const double& PRLPerilunePseudoLat,
		const double& PRLPerilunePseudoLong, const double& PRLPeriluneEccentric, const double& PRLPeriluneVelInc,
		CCoord3& posMJ2000, CCoord3& velMJ2000);
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
	/// @Param	LSOModElemEJ2000      ��������Ӱ����ʱ��EJ2000�������������
	/// @Return		
	//***********************************************************
	void CalInitPseudoPRL2LSO(const CCTime& timePRL, const double& PRLHeight, const double& PRLPerilunePseudoLat, const double& PRLPerilunePseudoLong,
		const double& PRLPeriluneEccentric, double PRLPeriluneVelInc, CModOrbElem& LSOModElemEJ2000);
	void CalInitElementsPRL2LSO(const CCTime& timePRL, const double& PRLHeight, const double& PRLPeriluneEccentric, const double& PRLPeriluneInc,
		const double& PRLPeriluneRAAN, double PRLPeriluneLAOP, CModOrbElem& LSOModElemEJ2000);
	//***********************************************************
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


