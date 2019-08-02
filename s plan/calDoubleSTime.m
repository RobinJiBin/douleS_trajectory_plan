function [ t ] = calDoubleSTime( v_s, v_e ,F ,D_max ,A_max ,J ,S )
%UNTITLED6 ����S���߹滮�ĸ��ε�ʱ�䣬
% ʱ���Ϊ T1,T2,T3,T4,T5,T6,T7
%   ����ֱ�Ϊ��
%        v_s    ��ʼ�ٶ�
%        v_e    �����ٶ�
%        F      �������ٶ�
%        D_max   �����ٶ�
%        A_max   �����ٶ�
%        J      ���Ӽ��ٶ�
%        S      ���о���

% ��һ������T1 - T6 ��
if F >= v_s + A_max^2/J
    T1 = A_max/J;
    T2 = (F - v_s)/A_max - A_max/J;
else
    T1 = sqrt((F - v_s)/J);
    T2 = 0;
end

if F >= v_e + D_max^2/J
    T5 = D_max/J;
    T6 = (F - v_e)/D_max - D_max/J;
else
    T5 = sqrt((F - v_e)/J);
    T6 = 0;
end
T3 = T1;
T7 = T5;

v1 = v_s + 0.5*J*T1^2;
v2 = v1 + J*T1*T2;
v3 = v_s + J*T1^2 + J*T1*T2;
v4 = v3;
v5 = v4 - 0.5*J*T5^2;
v6 = v5 - J*T5*T6;

S1 = (F + v_s)/2*(2*T1 + T2) + (F + v_e)/2*(2*T5 + T6);

if S1 < S
    T4 = (S - S1)/F; 
elseif S1 == S
    T4 = 0;
else
%     S1 > S ˵��û�����ٶΣ�T4 = 0������ϵͳû�д��������Ľ����ٶ�F
%     ���¼�������׶ε�����ʱ��
    T4 = 0;
    if v_s < v_e
%         ˵�����ٽ׶δ��ڼ��ٽ׶�ʱ��
%         ���տ��Ե�������ٶ�F1�����¼���
        F1 = v_e + D_max^2/J;
        F = F1;
        S1 = (F + v_s)/2*(2*T1 + T2) + (F + v_e)/2*(2*T5 + T6);
        if S1 < S
%             ˵��ϵͳ���ȼ���ʱ�䣬ʵ��S��Ϊ6��
              F = -A_max^2/2/J + sqrt(A_max^4 - 2*J*( A_max^2*(v_s + v_e) - J*(v_s^2 + v_e^2) - 2*A_max*J*S))/2/J;
              T1 = A_max/J;
              T3 = T1;
              T2 = (F - v_s )/A_max - A_max/J;
              T5 = D_max/J;
              T7 = T5;
              T6 = (F - v_e)/D_max - D_max/J;
        elseif S1 == S
%             ˵�������ȼ��ٶΣ�F1 == F
              T1 = A_max/J;
              T3 = T1;
              T2 = (F - v_s )/A_max - A_max/J;
              T5 = D_max/J;
              T7 = T5;
              T6 = 0;
        else 
%           S1 > S ˵�������ȼ��ٶȣ�ʵ������ٶ�С��F1        
%             ���F2��ȷ���Ƿ�����ȼ��ٶȶ�
            F2 = v_s + A_max^2/J;
            if F2 <= v_e
%                 �����ȼ��ٶΣ�S����Ϊ5��
%                 ���ݹ�ʽ11��������ȡ����ٶ�F
                F = calEquation11(v_s, v_e ,F ,D_max ,A_max ,J ,S);
                T1 = A_max/J;
                T3 = T1;
                T2 = (F - v_s )/A_max - A_max/J;
                T5 = sqrt(abs(F - v_e)/J);
                T7 = T5;
                T6 = 0;
            else
%                 F2 > ve �� S2��ȷ���Ƿ�����ȼ��ٶ�
                  F = F2;
                  T1 = A_max/J;
                  T5 = D_max/J;
                  S2 = (v_s + F)*T1 + (v_e + F)*T5;
                  if S2 >= S
%                       �������ȼ��ٶ�
%                         ʹ�÷���ʽ14�������F
                  F = calEquation14(v_s, v_e ,F ,D_max ,A_max ,J ,S);
                  T1 = sqrt((F - v_s)/J);
                  T3 = T1;
                  T2 = 0;
                  T6 = 0;
                  T5 = sqrt((F - v_e)/J);
                  T7 = T5;
                  else
%                       ���ݷ���11 �������� F��
                    F = calEquation11(v_s, v_e ,F ,D_max ,A_max ,J ,S);
                    T1 = A_max/J;
                    T3 = T1;
                    T2 = (F - v_s )/A_max - A_max/J;
                    T5 = sqrt(abs(F - v_e)/J);
                    T7 = T5;
                    T6 = 0;
                  end
            end
        end
    else
%     vs >= ve 
%       ��F1,����F1���S1
      F1 = v_s + A_max^2/J;
      F = F1;
      S1 = (F + v_s)/2*(2*T1 + T2) + (F + v_e)/2*(2*T5 + T6);
      if S1 < S
%           �����ȼ��ٶΣ����ڼ���ʱ��С�ڼ���ʱ�䣬��ض������ȼ��ٶ�
          F = -A_max^2/2/J + sqrt(A_max^4 - 2*J*( A_max^2*(v_s + v_e) - J*(v_s^2 + v_e^2) - 2*A_max*J*S))/2/J;
          T1 = A_max/J;
          T3 = T1;
          T2 = (F - v_s )/A_max - A_max/J;
          T5 = D_max/J;
          T7 = T5;
          T6 = (F - v_e)/D_max - D_max/J;        
      elseif S1 == S
          T1 = A_max/J;
          T3 = T1;
          T2 = 0;
          T5 = D_max/J;
          T7 = T5;
          T6 = (F - ve)/D_max - D_max/J;          
      else
%        S1 > S ˵���������ȼ��ٶΣ��������ж��Ƿ����ȼ��ٵ�
          F2 = v_e + D_max^2/J;
          F = F2;
          T1 = A_max/J;
          T5 = D_max/J;
          S2 = (v_s + F)*T1 + (v_e + F)*T5;
          if S2 >= S
%               �������ȼ��ٶΣ��������ȼ��ٶΣ�
%               ���ݷ���16�������ٶ�F������16�뷽��14��ͬ
                F = calEquation14(v_s, v_e ,F ,D_max ,A_max ,J ,S);
                T1 = sqrt((F - v_s)/J);
                T3 = T1;
                T2 = 0;
                T6 = 0;
                T5 = sqrt((F - v_e)/J);
                T7 = T5;
          else
%               ϵͳ�����ȼ��ٶΣ����ǲ������ȼ��ٶ�
%               ���ݷ���18�������ٶ�F
                 F = calEquation18(v_s, v_e ,F ,D_max ,A_max ,J ,S);
                 T1 = sqrt((F - v_s)/J);
                 T3 = T1;
                 T2 = 0;
                 T5 = D_max/J;
                 T7 = T5;
                 T6 = (F - v_e)/D_max - D_max/J;
          end
      end
    end     
end

    t(1) = T1;
    t(2) = T2;
    t(3) = T3;
    t(4) = T4;
    t(5) = T5;
    t(6) = T6;
    t(7) = T7;


end

