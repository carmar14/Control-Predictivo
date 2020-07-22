clc
close all
clear


%Matrices A,B,Q
Ts=0.1
Ad=[4/3 -2/3; 1 0];
Bd=[1 0]';
% Q=[4/9+0.001 -2/3; -2/3 1.001];
Q=[1 0; 0 1];
C=[-2/3 1];
R=1;
x0=[-0.5 1.3]'; %Condiciones inicales
ubounds=1; % Valores absolutos de las restricciones de u
xbounds=2;


% Ts=0.1; %sampling time
% Ad=[1.7 -0.4;0.4 0.8]; % discrete time system xk+1=Axk+Buk
% Bd=[0.2;0.05];          %                      yk=Cxk+Duk
% Cd=eye(2);
% Dd=0;
% Q=diag([1;1]); % Q and R matrices of sum x'Qx+u'Ru
% R=1; 
% x0=[0.75;-0.33]


% initial condition
% x0=[0.65;-0.33];%[-1,-1];%[0.07501;-0.3327];%[0.3;0.1];

%simulation length
simul_length=300;

% steady inputs
%%Us=[-0.8 -0.75 -0.65 -0.55 -0.45 -0.35];
%Us=[0 0 0 0 0 0];
% associated steady states
%%Xs=inv(eye(2)-Ad)*Bd*Us;
%Xs=inv(eye(2)-Ad)*Bd*Us;

% constraints
ubounds=1; % bound on the absolute value of u
xbounds=2;

% MPC parameters
horizon=5;
flag_with_preview=1;
% Q=diag([1;1]); % Q and R matrices of sum x'Qx+u'Ru
% R=1; % please change here to see the effect of different values of Q and R

% LQR solution used for the terminal penalty and terminal constraint
[K,P]=dlqr(Ad,Bd,Q,R); % terminal penalty for stability

% vectorial formulation
Qbig=kron(eye(horizon),Q);
Rbig=kron(eye(horizon),R);
%Qbig=blkdiag(Q,R);
Qbig=blkdiag(Qbig,Rbig);
% Rbig=kron(eye(horizon),R);
Abig=[];
Abig=kron(eye(horizon-1),Ad);
[m n]=size(Ad);
% Abig=[zeros(2,n+2); Abig];
A=zeros(horizon*m,horizon*m); 
A(m+1:end,1:end-n)=Abig;
Abig=A;

Bbig=[];
Bbig=kron(eye(horizon-1),Bd);
[m n]=size(Bd);
% Abig=[zeros(2,n+2); Abig];
B=zeros(horizon*m,horizon*n); 
B(m+1:end,1:end-n)=Bbig;
Bbig=B;

Hbig=[];
Hbig=eye(m);
H=zeros(horizon*m,m);
H(1:m,1:m)=Hbig;
Hbig=[H ]; %zeros(size(Bbig,2),size(Ad,2))];

% Bbig=0*kron(zeros(horizon+1,horizon+1),Bd);
% for i=0:horizon
%     Abig=[Abig;Ad^i];
% end
% for i=1:horizon+1
%     Bbig=Bbig+kron(diag(ones(horizon-i+1,1),-i),Ad^(i-1)*Bd);
% end
% Bbig=Bbig(:,1:end-size(Bd,2));

%cost function - we assume the presence of the quadratic term only
% H=Bbig'*Qbig*Bbig+Rbig;
H=Qbig;

%%Restricciones desigualdad
ExB=[Bbig;-Bbig];
A=kron(eye(horizon),[1;-1]);
Eu=A;
Ex=kron(eye(horizon),[1 0;0 -1]);
Ex=[Ex;Ex];
Ez=blkdiag(Ex,Eu);  %A para quadprog
bu=ubounds*ones(horizon,1);
xu=xbounds*ones(horizon*2,1);
Fz=blkdiag(xu,bu); % b para quadprog

%Restriccions de igualdad
[m n]=size(Abig);
I=eye(size(Abig,2));
Aeq=[I-Abig -Bbig];
beq=Hbig;



x1=[x0(1,1)];
x2=[x0(2,1)];
xx(1,:)=x0';
uu=[];
Tf=simul_length;
f=zeros(size(Ad,2)*horizon+size(Bd,2)*horizon,1);
f=f';
%f=zeros(3,horizon);
for ii=1:simul_length-1
    
   
    [U,FVAL,EXITFLAG]=quadprog(2*H,f,Ez,Fz);%,Aeq,beq*x0);
    if EXITFLAG==-2||isempty(U)
        flagtot=0;
    else
        x_t1=Ad*x0+Bd*U(horizon*size(Ad,2)+1); %Midiendo el nuevo estado a partir de la accion de control proveniente del problema de opt.
        x1=[x1 x_t1(1,1)];
        x2=[x2 x_t1(2,1)];
        x0=x_t1;
    end
end


% plotting

figure
t=0:Ts:Tf*Ts-Ts;
hold on
plot(t,x1,'-o')
plot(t,x2,'-ko')  
grid on
legend('x1 state','x2 state')




