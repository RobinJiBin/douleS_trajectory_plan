function [ value ] = equation18(  v_s, v_e ,F ,D_max ,A_max ,J ,S )
%UNTITLED2 ʹ��ţ������ѷ���������㷽��ʽ18�Ľ�
% �ú���Ϊ���̵����⣬���������ʽ��ֵ
%   ����ʽΪ �� 
%  ( F + v_e )/2 * ( D_max/J + ( F - v_e )/D_max ) + ( F + v_s
%   )*sqrt(( F - v_s )/J) = S
% ����FΪ��������
% v_s  v_e J A_max S Ϊ��֪����

value =  ( F + v_e )/2 * ( D_max/J + ( F - v_e )/D_max ) + ( F + v_s)*sqrt( abs(F - v_s )/J) - S;
end

