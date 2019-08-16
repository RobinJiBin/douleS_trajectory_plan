/********************************************************************************
* @file:  doubleSonline.cpp
* @note:  double S trajectory online calculation
* @auto:  jibin
* @other: 
********************************************************************************/
#include "doubleSonline.h"

S2TrajectoryPlan::S2TrajectoryPlan(): m_interval(0.02)
{
}

int S2TrajectoryPlan::GetNextMotionParam( stInstanceParam* TrajS, stInstanceParam* TrajE, stInstanceParam* curPara, stInstanceParam* prevPara )
{

	static double Td(0.0), Tj2a(0.0), Tj2b(0.0), hk(0.0);
	double t1(0.0), t2(0.0), t3(0.0), t4(0.0);

	//起始和终止点参数
	double q0 = TrajS->q;
	double q1 = TrajE->q;
	double v0 = TrajS->v;
	double v1 = TrajE->v;
	double a0 = TrajS->a;
	double a1 = TrajE->a;


	double vmax = m_MotionParams.vmax;
	double vmin = m_MotionParams.vmin;
	double amax = m_MotionParams.amax;
	double amin = m_MotionParams.amin;
	double jmax = m_MotionParams.jmax;//mm/s^3%j无穷大，速度梯形曲线
	double jmin = m_MotionParams.jmin;


	//k点--当前点;k_1点--前一个点
	double qk 	= curPara->q;
	double qk_1 = prevPara->q;
	double vk 	= curPara->v;
	double vk_1 = prevPara->v;
	double ak 	= curPara->a;
	double ak_1 = prevPara->a;
	double jk 	= curPara->j;
	double jk_1 = prevPara->j;
	

	

	const double Ts = m_interval;//单位s%时间间隔20ms
	static double Tk = 0;
	static double T_ = 0; 

	//phase1-加速和最大速；phase2-减速阶段
	//进入循环当，当前点没有接近终点时，循环继续

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
	//限幅
	if (ak>amax) ak = amax;
	else if (ak < amin)	ak = amin;
	
	if (jk<0&&ak_1>0&&ak<0)//ak在关键点的判断%避免震荡?
		ak = 0;
	//ac.push_back(ak);

	//根据ak值计算vk
	vk =  vk_1 + Ts*(ak + ak_1)/2;
	
    if (vk > vmax)	vk = vmax;
    if (vk < vmin) 	vk = vmin;//0
	//vc.push_back(vk);

	//根据vk值计算qk
	qk =  qk_1 + Ts*(vk + vk_1)/2; 
	//qc.push_back(qk);
	#endif

	//参数更新
	Tk = Tk + Ts;
	
	curPara->q = qk;
	curPara->v = vk;
	curPara->a = ak;
	curPara->j = jk;


	return 1;
}


int main(int argc, char* argv[])
{
	//输入q1、v0
	const int parain(2);
    int para[5]={0};
    if( (parain+1)!=argc )
    {
        cout<<parain<<" params are needed..."<<endl;
        cout<<"input: q1(mm), v0(mm/s)"<<endl;
        return -1;
    }
    else
    {
    	int i(0);
        for(i=0;i<parain;i++)
        	para[i]=atoi(argv[i+1]);
    }

    

	//定义记录曲线用变量
	vector<double> qc;
	vector<double> vc;
	vector<double> ac;
	vector<double> jc;


	S2TrajectoryPlan TraPlan;
	//当前点和上一时刻点参数
	stInstanceParam curPara, prevPara;
	//起始和终止点参数
	stInstanceParam TrajStart(0,0,0,0), TrajEnd(3500,0,0,0);
	TrajStart.v = para[0];//v0
	TrajEnd.q   = para[1];//q1


	int i(0);
	while(1)//(curPara.q <= q1)
	{
		int res = TraPlan.GetNextMotionParam( &TrajStart, &TrajEnd, &curPara, &prevPara);
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