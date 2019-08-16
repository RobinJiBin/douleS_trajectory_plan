/********************************************************************************
* @file:  doubleSonline.cpp
* @note:  double S trajectory online calculation
* @auto:  jibin
* @other: 
********************************************************************************/
#include "doubleSonline.h"
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
	double q;
	double v;
	double a;
	double j;
};

class S2TrajectoryPlan
{
public:
	S2TrajectoryPlan(){}
	virtual ~S2TrajectoryPlan(){}

public:
	int GetNextMotionParam( double curVel, double destDist, stInstanceParam* curPara, stInstanceParam* prevPara );



private:
	stMotionParam   m_MotionParams;
	stInstanceParam m_CurParam, m_prevParam;

	//起始和终止点参数
	double q0;
	double q1;//2000;%mm
	double v0;//800;  %%一直q1和v0求其他
	double v1;
	double a0;//mm/s^2
	double a1;


	//phase1-加速和最大速；phase2-减速阶段
	double phase_2 = 0; 
	double Tk = 0;//current time
	double Ts = 0.02;//单位s //时间间隔20ms
	double T_ = 0; 

};

int S2TrajectoryPlan::GetNextMotionParam( double curVel, double destDist, stInstanceParam* curPara, stInstanceParam* prevPara )
{

	static double Td(0.0), Tj2a(0.0), Tj2b(0.0), hk(0.0);
	double t1(0.0), t2(0.0), t3(0.0), t4(0.0);

	//起始和终止点参数
	double q0 = 0;
	double q1 = 3500;//2000;%mm
	double v0 = 0;//800;  %%一直q1和v0求其他
	double v1 = 0;
	double a0 = 0;//mm/s^2
	double a1 = 0;


	double vmax = 800;//mm/s
	double vmin = 0;
	double amax = 250;//mm/s^2
	double amin = -250;
	double jmax = 400;//mm/s^3%j无穷大，速度梯形曲线
	double jmin = -400;
	//double vmax_org = vmax;%??


	//k点--当前点;k_1点--前一个点
	double qk 	= curPara->q;
	double qk_1 = prevPara->q;
	double vk 	= curPara->v;
	double vk_1 = prevPara->v;
	double ak 	= curPara->a;
	double ak_1 = prevPara->a;
	double jk 	= curPara->j;
	double jk_1 = prevPara->j;
	

	//phase1-加速和最大速；phase2-减速阶段
	double phase_2 = 0; 
	static double Tk = 0;
	double Ts = 0.02;//单位s%时间间隔20ms
	static double T_ = 0; 

		//进入循环当，当前点没有接近终点时，循环继续
	//while (1)
	{
		cout<<"T_ = "<<T_<<endl;
	    //判断是否已经进入到减速阶段
	    if (T_>0) //通过T_判断是否进入减速阶段j
	    {
	        t1 = Tk - T_;
	        t2 = Tj2a;
	        t3 = (Td - Tj2b);
	        t4 = Td;
	        //已经进入到减速阶段
	        //直接判断jk
	        if ((Tk - T_)>=0&&(Tk - T_)<(Tj2a)){
	        	cout<<"1"<<endl;
	         	jk = jmin;
	         }            
	        else if( (Tk - T_)>=(Tj2a)&&(Tk - T_)<(Td - Tj2b) ){
	        	cout<<"2"<<endl;
	            jk = 0;
	        }             
	        else if((Tk - T_)>=((Td - Tj2b))&&(Tk - T_)<(Td))
	        {
	             cout<<"3"<<endl;
	             jk = jmax;
	        }            
	        else{        	
	            cout<<"4"<<endl;
	            //break; //退出位置
	            return 0;
	        }

	        //Tk - T_//disp?
	    }
	    else{
	        //未进入到减速阶段
	        //计算Td,Tj2a,Tj2b，用以计算hk。
	        Tj2a = (amin - ak)/jmin;
	        Tj2b = (a1 - amin)/jmax;
	        Td = (v1 - vk)/amin + Tj2a*(amin - ak)/2/amin + Tj2b*(amin - a1)/2/amin;//原理？

	        //首先判断在减速段，是否能够达到最小减速度amin
	        if (Td <= Tj2a + Tj2b){
	        	double tmp1 = ((jmax - jmin)*(ak*ak*jmax - jmin*(a1*a1 + 2*jmax*(vk - v1))));
	        	tmp1 = pow(tmp1,2);
	        	Tj2a = -1*ak/jmin + tmp1/jmin/(jmin - jmax);
	           	Tj2b = a1/jmax    + tmp1/jmin/(jmin - jmax);
	           	Td = Tj2a + Tj2b;
	        }           

	        //计算hk，判断是否需要进入减速阶段
	        hk = 0.5*ak*Td*Td + (jmin*Tj2a*(3*Td*Td - 3*Td*Tj2a + Tj2a*Tj2a) + jmax*pow(Tj2b,3))/6 + Td*vk;
	        if (hk <(q1 - qk)){
	            //case 1 加速或匀速运动段
	            //判断加加速度的值
	            if ( ((vk - ak*ak/2/jmin) < vmax)&&(ak < amax) ){
	            	//disp('##############4##############')                
	                jk = jmax;
	            }
	            else if( ((vk - ak*ak/2/jmin) < vmax)&&(ak >= amax) ){
	            	//disp('##############5##############')              
	                jk = 0;
	            }                 
	            else if( ((vk - ak*ak/2/jmin) >= vmax)&&(ak >0) ){
	            	//disp('##############6##############')                
	                jk = jmin;
	            }          
	            else if( ((vk - ak*ak/2/jmin) >= vmax)&&(ak <=0) ){
	            	//disp('##############7##############'
	                jk = 0;
	            }
	        }
	        //case 2 减速度阶段
	        else{
	        	//进入到减速阶段，记录时间
	            T_ = Tk;
	            jk = jmin;
	            //disp('8')
	        }
	            
	    }//end if T_>0
	    
	    #if 1
		//根据jk计算本周期的加速度，速度，位置

	   	//jc.push_back(jk);
	   
		ak = ak_1 + Ts*(jk + jk)/2 ; 
		if (ak>amax) ak = amax;
		else if (ak < amin)	ak = amin;
		
		if (jk<0&&ak_1>0&&ak<0)//ak在关键点的判断%避免震荡?
			ak = 0;

		//ac.push_back(ak);

		//根据ak值计算vk
		vk =  vk_1 + Ts*(ak + ak_1)/2;
		/*
		%    if vk> vmax
		%        vk = vmax;
		%    end
		%    if vk < vmin
		%        vk = vmin;
		%    end*/

		//vc.push_back(vk);

		//根据vk值计算qk
		qk =  qk_1 + Ts*(vk + vk_1)/2; 
		
		//qc.push_back(qk);
		#endif

		//参数更新
		qk_1 = qk;
		vk_1 = vk;
		ak_1 = ak;
		jk_1 = jk;
		Tk = Tk + Ts;

		
		curPara->q = qk;
		curPara->v = vk;
		curPara->a = ak;
		curPara->j = jk;


		return 1;
	   
	}  //end while



}


int main(int argc, char* argv[])
{
	//输入v0、q1
	const int parain(2);
    int para[5]={0};
    if( (parain+1)!=argc )
    {
        cout<<parain<<" params are needed..."<<endl;
        cout<<" v0(mm/s); q1(mm)"<<endl;
        return -1;
    }
    else
    {
    	int i(0);
        for(i=0;i<parain;i++)
        	para[i]=atoi(argv[i+1]);
    }

    

    //起始和终止点参数
	double q0 = 0;
	double q1 = 3500;//para[1];//2000;%mm
	double v0 = 0;//para[0];//800;  %%一直q1和v0求其他
	double v1 = 0;
	double a0 = 0;//mm/s^2
	double a1 = 0;
	//定义记录曲线用变量
	vector<double> qc;
	vector<double> vc;
	vector<double> ac;
	vector<double> jc;


	S2TrajectoryPlan TraPlan;
	stInstanceParam curPara, prevPara;


	int i(0);
	while(1)//(curPara.q <= q1)
	{
		int res = TraPlan.GetNextMotionParam( v0, q1, &curPara, &prevPara);
		if(!res) break;//减速阶段完成break；

		vc.push_back(curPara.v);
		qc.push_back(curPara.q);
		ac.push_back(curPara.a);
		jc.push_back(curPara.j);
		
		cout<<"prevPara.[q.v.a.j] = "<<prevPara.q<<", "<<prevPara.v<<", "<<prevPara.a<<", "<<prevPara.j<<endl;
		cout<<"curPara.[q.v.a.j] = "<<curPara.q<<", "<<curPara.v<<", "<<curPara.a<<", "<<curPara.j<<endl;
	
		prevPara = curPara;
		i++;
		cout<<"i = "<<i<<endl;
		cout<<"q1 = "<<q1<<endl;
		//if(i>220) break;
	}



	ofstream outFile;
    outFile.open("qc.txt");
    outFile<<fixed;
    outFile.precision(3);
    outFile.setf(ios_base::showpoint);
	vector<double>::iterator it;
	for(it=qc.begin(); it!=qc.end(); it++)
		outFile<<*it<<endl;
	cout<<"qc.size: "<<qc.size()<<endl;
	outFile.close();

	outFile.open("vc.txt");
	for(it=vc.begin(); it!=vc.end(); it++)
		outFile<<*it<<endl;
	cout<<"vc.size: "<<vc.size()<<endl;
	outFile.close();

	outFile.open("ac.txt");
	for(it=ac.begin(); it!=ac.end(); it++)
		outFile<<*it<<endl;
	cout<<"ac.size: "<<ac.size()<<endl;
	outFile.close();

	outFile.open("jc.txt");
	for(it=jc.begin(); it!=jc.end(); it++)
		outFile<<*it<<endl;
	cout<<"jc.size: "<<jc.size()<<endl;
	outFile.close();


	return 0;
}