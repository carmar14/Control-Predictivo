clear
clc
close all

%AUTOR: CARLOS MARIO PAREDES
%Matrices A,B,Q
A=[4/3 -2/3; 1 0];
B=[1 0]';
Q=[4/9+0.001 -2/3; -2/3 1.001];
C=[-2/3 1];
R=0.001;
Pf=Q;
N=6;   %5+1 por la programacion de matlab
PI=zeros(2,2,N);
PI(:,:,end)=Pf;
uo=zeros(1,N-1);   %Entrada-accion de control optima


%Calculo de PI
for k=N:-1:2
    PI(:,:,k-1)=Q+A'*PI(:,:,k)*A-A'*PI(:,:,k)*B*inv((B'*PI(:,:,k)*B+R))*B'*PI(:,:,k)*A;
end

%Calculo de la  ganacia K
K=zeros(1,2,N-1);
for k=N:-1:2
    K(:,:,k-1)=-inv(B'*PI(:,:,k)*B+R)*B'*PI(:,:,k)*A;
end
x0=[0.5 0.5]';%Condiciones iniciales
xk=zeros(2,1,N-1);
xk(:,:,1)=x0;
%Calculo de accion de control optimo
for k=1:N-1
    xk(:,:,k+1)=A*xk(:,:,k)+B*K(:,:,k)*xk(:,:,k);
    u(k)=K(:,:,k)*xk(:,:,k);
end

%Costo optimo
Vk=zeros(1,N);
for k=1:N
    Vk(k)=0.5*xk(k)*PI(k)*xk(k);
end


%Autovalores en lazo cerrado
eig(A+B*K(:,:,1))  %Autovalores para N=5, sistema inestable

[Kc,S,e] = dlqr(A,B,Q,R,[0 0]') 
eig(A-B*Kc) %Autovalores con ganancia calculada por dlqr