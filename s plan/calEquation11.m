function [ F ] = calEquation11( v_s, v_e ,F ,D_max ,A_max ,J ,S)
%UNTITLED2 ����ţ������ѷ���������㷽��ʽ11��ֵ
%   
F0 = 5;
for i = 1:1:10
    v0 =  equation11(v_s, v_e ,F0 ,D_max ,A_max ,J ,S);
    dot0 = dotEquation11(v_s, v_e ,F0 ,D_max ,A_max ,J ,S);
    F1 = F0 - v0/dot0;
    v1 = equation11(v_s, v_e ,F1 ,D_max ,A_max ,J ,S);
    
    if abs(F1 - F0) < 0.001|| abs(F1) < 0.001
        disp('��������');
        F1;
        break;
    end
    F0 = F1;
end
F = F1;
end

