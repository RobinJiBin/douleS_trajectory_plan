function [ value ] = dotEquation14(  v_s, v_e ,F ,D_max ,A_max ,J ,S)
%UNTITLED2 ʹ��ţ������ѷ���������㷽��ʽ11�Ľ�
% �ú���Ϊ���̵����⣬���������ʽ��ֵ
%   ����ʽΪ �� 
%  ( F + v_s )*sqrt((F - v_s)/J) + (F + v_e)*sqrt((F - v_e)/J) - S
%  ���̵ĵ���Ϊ��
%  ((F - v_e)/J)^(1/2) + ((F - v_s)/J)^(1/2) + (F + v_e)/(2*J*((F - v_e)/J)^(1/2)) + (F + v_s)/(2*J*((F - v_s)/J)^(1/2))
% ����FΪ��������
% v_s  v_e J A_max S Ϊ��֪����

value = (abs(F - v_e)/J)^(1/2) + (abs(F - v_s)/J)^(1/2) + (F + v_e)/(2*J*(abs(F - v_e)/J)^(1/2)) + (F + v_s)/(2*J*(abs(F - v_s)/J)^(1/2));
end

