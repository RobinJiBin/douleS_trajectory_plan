function [ value ] = equation14(  v_s, v_e ,F ,D_max ,A_max ,J ,S )
%UNTITLED2 ʹ��ţ������ѷ���������㷽��ʽ14�Ľ�
% �ú���Ϊ���̵����⣬���������ʽ��ֵ
%   ����ʽΪ �� 
%  ( F + v_s )*sqrt((F - v_s)/J) + (F + v_e)*sqrt((F - v_e)/J) - S
%   
% ����FΪ��������
% v_s  v_e J A_max S Ϊ��֪����

value = ( F + v_s )*sqrt(abs(F - v_s)/J) + (F + v_e)*sqrt(abs(F - v_e)/J) - S;
end

