%��д����ģ�����߼���˫S���߼Ӽ���
%�ӳ�ʼ��q0��ʼ��һֱ���㵽�յ�q1����
%�����Ҫ�ı���

q0 = 0;
q1 = 10;
v0 = 0;
v1 = 3;
a0 = 1;
a1 = 0;


vmax = 5;
vmin = -5;
amax = 10;
amin = -8;
jmax = 30;
jmin = -40;
vmax_org = vmax;
%������k���ϵ�
qk = q0;
qk_1 = q0;
vk = v0;
vk_1 = v0;
ak = a0;
ak_1 = a0;
jk = 0;
jk_1 = 0;
e_min = 0.001;


phase_2 = 0; 
Tk = 0;
Ts = 0.001;
T_ = 0; 

%�����¼�����ñ���
qc = [];
vc = [];
ac = [];
jc = [];

%����ѭ��������ǰ��û�нӽ��յ�ʱ��ѭ������
while 1
%     if qk > 2.5 && qk < 4
%         vmax = vmax_org*0.8;
%     end
    %�ж��Ƿ��Ѿ����뵽���ٽ׶�

    if T_>0
        t1 = Tk - T_;
        t2 = Tj2a;
        t3 = (Td - Tj2b);
        t4 = Td;
        %�Ѿ����뵽���ٽ׶�
    %ֱ���ж�jk
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
        %δ���뵽���ٽ׶�
        %����Td,Tj2a,Tj2b�����Լ���hk��
        Tj2a = (amin - ak)/jmin;
        Tj2b = (a1 - amin)/jmax;
        Td = (v1 - vk)/amin + Tj2a*(amin - ak)/2/amin + Tj2b*(amin - a1)/2/amin;

        %�����ж��ڼ��ٶΣ��Ƿ��ܹ��ﵽ��С���ٶ�amin
        if Td <= Tj2a + Tj2b
           Tj2a = -1*ak/jmin + ((jmax - jmin)*(ak^2*jmax - jmin*(a1^2 + 2*jmax*(vk - v1))))^2/jmin/(jmin - jmax);
           Tj2b = a1/jmax + ((jmax - jmin)*(ak^2*jmax - jmin*(a1^2 + 2*jmax*(vk - v1))))^2/jmin/(jmin - jmax);
           Td = Tj2a + Tj2b;
        end
         %����hk���ж��Ƿ���Ҫ������ٽ׶�
        hk = 0.5*ak*Td^2 + (jmin*Tj2a*(3*Td^2 - 3*Td*Tj2a + Tj2a^2) + jmax*Tj2b^3)/6 + Td*vk;
        if hk <(q1 - qk)
            %case 1 ���ٻ������˶���
            %�жϼӼ��ٶȵ�ֵ
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
        %case 2 ���ٶȽ׶�
        else
            %���뵽���ٽ׶Σ���¼ʱ��
            T_ = Tk;
            jk = jmin;
%             disp('8')
        end
    end%end if T_>0
    
%����jk���㱾���ڵļ��ٶȣ��ٶȣ�λ��

   jc = [jc jk];
   
   ak = ak_1 + Ts*(jk + jk)/2 ; 
   if ak>amax
       ak = amax;
   end
   if ak < amin
       ak = amin;
   end
   if jk<0&&ak_1>0&&ak<0
       ak = 0;
   end  
   ac = [ac ak]; 
   vk =  vk_1 + Ts*(ak + ak_1)/2;
%    if vk> vmax
%        vk = vmax;
%    end
%    if vk < vmin
%        vk = vmin;
%    end

       
   vc = [vc vk];
   qk =  qk_1 + Ts*(vk + vk_1)/2; 
    qc = [qc qk];

   qk_1 = qk;
   vk_1 = vk;
   ak_1 = ak;
   jk_1 = jk;
   Tk = Tk + Ts;
   
end%end while

plot(qc)
figure
plot(vc)
figure
plot(ac)
figure
plot(jc)
