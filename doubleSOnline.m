%编写程序模拟在线计算双S曲线加减速
%从初始点q0开始，一直计算到终点q1结束
%定义必要的变量

%起始和终止点参数
q0 = 0000;
q1 = 3500;%mm
v0 = 0;  %%一直q1和v0求其他
v1 = 0;
a0 = 0;%mm/s^2
a1 = 0;


vmax = 800;%mm/s
vmin = 0;
amax = 250;%mm/s^2
amin = -250;
jmax = 400;%mm/s^3%j无穷大，速度梯形曲线
jmin = -400;
vmax_org = vmax;%??


%定义在k点上的%k点是哪个点
qk = q0;
qk_1 = q0;
vk = v0;
vk_1 = v0;
ak = a0;
ak_1 = a0;
jk = 0;
jk_1 = 0;
%e_min = 0.001;%?


phase_2 = 0; 
Tk = 0;
Ts = 0.02%时间间隔20ms
T_ = 0; 

%定义记录曲线用变量
qc = [];
vc = [];
ac = [];
jc = [];

%进入循环当，当前点没有接近终点时，循环继续
while 1
%     if qk > 2.5 && qk < 4
%         vmax = vmax_org*0.8;
%     end
    %判断是否已经进入到减速阶段

    if T_>0 %通过T_判断是否进入减速阶段?
        t1 = Tk - T_;
        t2 = Tj2a;
        t3 = (Td - Tj2b);
        t4 = Td;
        %已经进入到减速阶段
    %直接判断jk
        if (Tk - T_)>=0&&(Tk - T_)<(Tj2a)
             disp('1')
            jk = jmin;
        elseif (Tk - T_)>=(Tj2a)&&(Tk - T_)<(Td - Tj2b)
             disp('2')
            jk = 0;
        elseif(Tk - T_)>=((Td - Tj2b))&&(Tk - T_)<(Td)
             disp('3')
            jk = jmax;
        else
            disp('4')
%             jk =0;
%             ak =a1;
%             vk =v1;
%             qk = q1;
            break;
        end
        Tk - T_
    else
        %未进入到减速阶段
        %计算Td,Tj2a,Tj2b，用以计算hk。
        Tj2a = (amin - ak)/jmin;
        Tj2b = (a1 - amin)/jmax;
        Td = (v1 - vk)/amin + Tj2a*(amin - ak)/2/amin + Tj2b*(amin - a1)/2/amin;

        %首先判断在减速段，是否能够达到最小减速度amin
        %%如果条件满足，重新计算Td
        if Td <= Tj2a + Tj2b
           Tj2a = -1*ak/jmin + ((jmax - jmin)*(ak^2*jmax - jmin*(a1^2 + 2*jmax*(vk - v1))))^2/jmin/(jmin - jmax);
           Tj2b = a1/jmax + ((jmax - jmin)*(ak^2*jmax - jmin*(a1^2 + 2*jmax*(vk - v1))))^2/jmin/(jmin - jmax);
           Td = Tj2a + Tj2b;
        end
        
         %计算hk，判断是否需要进入减速阶段
        hk = 0.5*ak*Td^2 + (jmin*Tj2a*(3*Td^2 - 3*Td*Tj2a + Tj2a^2) + jmax*Tj2b^3)/6 + Td*vk;
        if hk <(q1 - qk)
            %case 1 加速或匀速运动段
            %判断加加速度的值
            if ((vk - ak^2/2/jmin) < vmax)&&(ak < amax)
%                 disp('##############4##############')
                
                jk = jmax;
            elseif((vk - ak^2/2/jmin) < vmax)&&(ak >= amax)
%                 disp('##############5##############')
              
                jk = 0;
            elseif((vk - ak^2/2/jmin) >= vmax)&&(ak >0 )
%                 disp('##############6##############')
                
                jk = jmin;
            elseif((vk - ak^2/2/jmin) >= vmax)&&(ak <=0 )
                disp('##############7##############')
                vk
                ak
                jk = 0;
            end
        %case 2 减速度阶段
        else
            %进入到减速阶段，记录时间
            T_ = Tk;
            jk = jmin;
%             disp('8')
        end
    end%end if T_>0
    
%根据jk计算本周期的加速度，速度，位置

   jc = [jc jk];
   
   ak = ak_1 + Ts*(jk + jk)/2 ; 
   if ak>amax
       ak = amax;
   end
   if ak < amin
       ak = amin;
   end
   if jk<0&&ak_1>0&&ak<0%ak在关键点的判断%避免震荡?
       ak = 0;
   end  
   ac = [ac ak]; 
   
%根据ak值计算vk
   vk =  vk_1 + Ts*(ak + ak_1)/2;
%    if vk> vmax
%        vk = vmax;
%    end
%    if vk < vmin
%        vk = vmin;
%    end

   vc = [vc vk];

%根据vk值计算qk
   qk =  qk_1 + Ts*(vk + vk_1)/2; 
   qc = [qc qk];

   qk_1 = qk;
   vk_1 = vk;
   ak_1 = ak;
   jk_1 = jk;
   Tk = Tk + Ts;
   
end%end while

subplot(411)
plot(qc,'LineWidth',2)
ylabel('Q距离(mm)')%y轴标记
title('trajectory plan')%标题
grid
subplot(412)
plot(vc,'LineWidth',2)
grid
ylabel('V速度(mm/s)')%y轴标记
subplot(413)
plot(ac,'LineWidth',2)
grid
ylabel('A加速度(mm/s^2)')%y轴标记
subplot(414)
plot(jc,'LineWidth',2)
grid
ylabel('J加加速度(mm/s^3)','FontSize',12)%y轴标记
xlabel('时间(20ms)')%x轴标记








