function [ value ] = equation11(  v_s, v_e ,F ,D_max ,A_max ,J ,S )
%UNTITLED2 ʹ��ţ������ѷ���������㷽��ʽ11�Ľ�
% �ú���Ϊ���̵����⣬���������ʽ��ֵ
%   ����ʽΪ �� 
%  ( F + v_s )/2 * ( A_max/J + ( F - v_s )/A_max ) + ( F + v_e
%   )*sqrt(( F - v_e )/J) = S
% ����FΪ��������
% v_s  v_e J A_max S Ϊ��֪����

value = ( F + v_s )/2 * ( A_max/J + ( F - v_s )/A_max ) + ( F + v_e)*sqrt(abs( F - v_e )/J) - S;
end

