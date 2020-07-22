clc
close all
clear all

Ts=0.1; %sampling time
Ad=[4/3 -2/3;1 0]; % discrete time system xk+1=Axk+Buk
Bd=[1 0]';         %                      yk=Cxk+Duk

x0=[0.2;-0.1]; % initial condition of the system

ubounds=1; % bound on the absolute value of u
xbounds=2;

Q=([4/9+0.001 -2/3;-2/3 1.001]); % Q and R matrices of sum x'Qx+u'Ru
R=100; % please change here to see the effect of different values of Q and R

horizon=30;


%%


% below we construct all the matrices according to what described in the
% scripts

Qbig=kron(eye(horizon),Q);
Rbig=kron(eye(horizon),R);

Abig=[];
Bbig=0*kron(eye(horizon),Bd);

for i=1:horizon
    Abig=[Abig;Ad^(i-1)];
    Bbig=Bbig+kron(diag(ones(horizon-i,1),-i),Ad^(i-1)*Bd);
end



%cost function - we assume the presence of the quadratic term only

H=Bbig'*Qbig*Bbig+Rbig;

% constraints on the inputs - no state constraints for the moment. Please
% add them

A=kron(eye(horizon),[1;-1]);
b=ubounds*ones(horizon*2,1);

A=[A;Bbig;-Bbig];
b=[b;-Abig*x0+xbounds*ones(horizon*2,1);Abig*x0+xbounds*ones(horizon*2,1)];
f=Bbig'*Qbig*Abig*x0;


[U,FVAL,EXITFLAG]=quadprog(2*H,2*f,A,b);

if EXITFLAG==-2||isempty(U)
    flagtot=0;
else
    X=Abig*x0+Bbig*U;
    X_re=[reshape(X,2,horizon)'];
    
    figure(1)
    
    hold on
    subplot(1,2,1);plot(0:Ts:Ts*(horizon-1),X_re(:,1),'b');hold on;plot(0:Ts:Ts*(horizon-1),X_re(:,1),'bo');plot(0:Ts:Ts*(horizon-1),X_re(:,2),'r');plot(0:Ts:Ts*(horizon-1),X_re(:,2),'ro');grid on
    subplot(1,2,2);plot(0:Ts:Ts*(horizon-1),U,'k');hold on;plot(0:Ts:Ts*(horizon-1),U,'ko');grid on
    
    
end

