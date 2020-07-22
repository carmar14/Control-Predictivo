clear
close all
clc

%Codigo prueba
A=[0.9 1; 0 0.9];
B=[0.5; 1];
C=[1 0];

tot=[eye(2)-A -B; C 0];
it=inv(tot);
%b=[0 0 r]';
%[x1 x2 uss]'=it*b;
%x1=r;
%x2=0.0952;
%uss=0.0095;
