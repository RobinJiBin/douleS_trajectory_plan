/********************************************************************************
* @file:  doubleSonline.cpp
* @note:  double S trajectory online calculation
* @auto:  jibin
* @other: 
********************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <math.h>

using namespace std;


struct stMotionParam
{
	stMotionParam(): vmax(800),vmin(0),amax(250),amin(-250),jmax(400),jmin(-400) {}
	double vmax;//mm/s
	double vmin;
	double amax;//mm/s^2
	double amin;
	double jmax;//mm/s^3%j无穷大，速度梯形曲线
	double jmin;
};

struct stInstanceParam
{
	//k点--当前点;k_1点--前一个点
	stInstanceParam(): q(0),v(0),a(0),j(0) {}
	stInstanceParam(double q, double v, double a, double j): 
					q(0),v(0),a(0),j(0) {}
	double q;//mm
	double v;//mm/s
	double a;//mm/s^2
	double j;//mm/s^3
};

class S2TrajectoryPlan
{
public:
	S2TrajectoryPlan();
	virtual ~S2TrajectoryPlan(){}

public:
	int GetNextMotionParam( stInstanceParam* TrajS, stInstanceParam* TrajE, stInstanceParam* curPara, stInstanceParam* prevPara );



private:
	stMotionParam   m_MotionParams;
	stInstanceParam m_CurParam, m_prevParam;

	double m_interval;

};