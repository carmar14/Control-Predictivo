clear 
close all
clc

%Matrices A,B,Q
Ts=0.1
Ad=[4/3 -2/3; 1 0];
Bd=[1 0]';
Q=[4/9+0.001 -2/3; -2/3 1.001];
C=[-2/3 1];
R=0.001;
x0=[-0.5 1.3]'; %Condiciones inicales
ubounds=1; % Valores absolutos de las restricciones de u
xbounds=2;

%N=30; %Horizonte a 30 muestras
N=15;

Qbig=kron(eye(N),Q);
Rbig=kron(eye(N),R); 

Abig=[];
Bbig=0*kron(eye(N),Bd);
for i=1:N
    Abig=[Abig;Ad^(i-1)];
    Bbig=Bbig+kron(diag(ones(N-i,1),-i),Ad^(i-1)*Bd);
end

%funcion de costo
H=Bbig'*Qbig*Bbig+Rbig;

%Restricciones: Ax<=b para quadprog
%-1<=u<=1
%-1<=x1<=1    y     -1<=x2<=1 =====>  x1<=1, -x1<=1, -x2<=1, x2<=1;

% Control 
Tf=50;   %Extendiendo a 50 muestras
x1=[x0(1,1)];
x2=[x0(2,1)];
A=kron(eye(N),[1;-1]);
A=[A;Bbig;-Bbig];
for i=1:Tf

    b=ubounds*ones(N*2,1);
    b=[b;-Abig*x0+xbounds*ones(N*2,1);Abig*x0+xbounds*ones(N*2,1)];
    f=Bbig'*Qbig*Abig*x0;
    [U,FVAL,EXITFLAG]=quadprog(2*H,2*f,A,b);
    if EXITFLAG==-2||isempty(U)
        flagtot=0;
    else
        x_t1=Ad*x0+Bd*U(1); %Midiendo el nuevo estado a partir de la accion de control proveniente del problema de opt.
        x1=[x1 x_t1(1,1)];
        x2=[x2 x_t1(2,1)];
        x0=x_t1;
    end
    
    clc

end

figure
t=0:Ts:Tf*Ts;%-Ts;
hold on
plot(t,x1,'-o')
plot(t,x2,'-ko')  
grid on
legend('x1 state','x2 state')
    

