function [ value ] = dotEquation18(v_s, v_e ,F ,D_max ,A_max ,J ,S)
%UNTITLED2 ʹ��ţ������ѷ���������㷽��ʽ11�Ľ�
% �ú���Ϊ���̵����⣬���������ʽ��ֵ
%   ����ʽΪ �� 
%  ( F + v_s )/2 * ( A_max/J + ( F - v_s )/A_max ) + ( F + v_e
%   )*sqrt(( F - v_e )/J) = S
%  ���̵ĵ���Ϊ��
%  A_max/(2*J) + (F - v_s)/(2*A_max) + (F/2 + v_s/2)/A_max + ((F - v_e)/J)^(1/2) + (F + v_e)/(2*J*((F - v_e)/J)^(1/2))
%  D_max/(2*J) + (F - v_e)/(2*D_max) + (F/2 + v_e/2)/D_max + ((F - v_s)/J)^(1/2) + (F + v_s)/(2*J*((F - v_s)/J)^(1/2))
% ����FΪ��������
% v_s  v_e J A_max S Ϊ��֪����

value = D_max/(2*J) + (F - v_e)/(2*D_max) + (F/2 + v_e/2)/D_max + (abs(F - v_s)/J)^(1/2) + (F + v_s)/(2*J*(abs(F - v_s)/J)^(1/2));
end

