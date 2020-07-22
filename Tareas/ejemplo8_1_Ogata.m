clc
clear
close
%Ejemplo de control optimo libro de ogata Sistemas de control en tiempo
%discreto
A=0.3679;
B=0.6321;
x0=1;
N=9;
Q=1;  
R=1;
Pf=1;%S
PI=zeros(1,N);
PI(end)=Pf;
%Calculo de PI
for k=N:-1:2
    PI(k-1)=Q+A'*PI(k)*A-A'*PI(k)*B*inv((B'*PI(k)*B+R))*B'*PI(k)*A;
end

%Calculo de la  ganacia K
K=zeros(1,N);
for k=N:-1:2
    K(k-1)=-inv(B'*PI(k)*B+R)*B'*PI(k)*A;
end

%Calculo de accion de control optimo
u=zeros(1,N);
%Evolucion del estado
xk=zeros(1,N)';
xk(1)=x0;
%x_k1=A*xk+B*K*xk;
for k=1:N
    xk(k+1)=A*xk(k)+B*K(k)*xk(k);
    u(k)=K(k)*xk(k);
end


