clc
close all
clear


Ts=0.1; %sampling time
Ad=[1.7 -0.4;0.4 0.8]; % discrete time system xk+1=Axk+Buk
Bd=[0.2;0.05];          %                      yk=Cxk+Duk
Cd=eye(2);
Dd=0;
W = Polyhedron([ eye(2) ; -eye(2) ],0.1*[ 1 1 1 1 ]');

% initial condition
x0=[0.07501;-0.3327];%[-1,-1];%[0.07501;-0.3327];%[0.3;0.1];

%simulation length
simul_length=300;


% constraints
ubounds=6; % bound on the absolute value of u
Uset=Polyhedron([ eye(1) ; -eye(1) ],ubounds*[ 1 1 ]');


xbounds=5;
Xset=Polyhedron([ eye(2) ; -eye(2) ],xbounds*[ 1 1 1 1]');

% MPC parameters
horizon=15;
flag_with_preview=1;
Q=diag([1;1]); % Q and R matrices of sum x'Qx+u'Ru
R=1; % please change here to see the effect of different values of Q and R

% LQR solution used for the terminal penalty and terminal constraint
[K,P]=dlqr(Ad,Bd,Q,R); % terminal penalty for stability

% approximation for the mRPI set
epsilon = 5*10^-5;
F  = epsilon_mRPI(Ad-Bd*K,W,epsilon);
[HH,KK]=F.doubleHK;%Hx<=K
mRPIset=Polyhedron('A',HH,'b',KK);
Ured=Uset-K*mRPIset;
HKtemp=Ured.H;
UredH=HKtemp(:,1:size(Bd,2));
UredK=HKtemp(:,end);

Xred=Xset-mRPIset;
HKtemp=Xred.H;
XredH=HKtemp(:,1:size(Ad,2));
XredK=HKtemp(:,end);

[tempHHKK HH KK set_X]=fsetX2(Ad,Bd,K,Ured,Xred,zeros(size(Ad,2),1),zeros(size(Bd,2),1));


plot(Xset)
hold on
plot(Xred,'Color','y')
plot(set_X{end},'Color','g')
plot(mRPIset,'Color','b')


% vectorial formulation
Qbig=kron(eye(horizon),Q);
Qbig=blkdiag(Qbig,P);
Rbig=kron(eye(horizon),R);
Abig=[];
Bbig=0*kron(zeros(horizon+1,horizon+1),Bd);
for i=0:horizon
    Abig=[Abig;Ad^i];
end
for i=1:horizon+1
    Bbig=Bbig+kron(diag(ones(horizon-i+1,1),-i),Ad^(i-1)*Bd);
end
Bbig=Bbig(:,1:end-size(Bd,2));

%cost function - we assume the presence of the quadratic term only
H=Bbig'*Qbig*Bbig+Rbig;


xx(1,:)=x0';
uu=[];


for ii=1:simul_length-1
    
    
    x0=xx(ii,:)';
    
    % definition of the bounds on the input
%     A=kron(eye(horizon),[1;-1]);
%     b=ubounds*ones(horizon*2,1);

    A=kron(eye(horizon),UredH);
    b=kron(ones(horizon,1),UredK);


    % temp and lastline are used to impose the terminal constraint 
    temp=Abig*x0;
    lastline=temp(end-1:end,:);
    
    % the two lines below include the previous input constraints, the state
    % constraints expressed in terms of inputs and the terminal constraint
%     A=[A;Bbig;-Bbig;HH*Bbig(end-1:end,:)];
%     b=[b;-Abig*x0+xbounds*ones((horizon+1)*2,1);Abig*x0+xbounds*ones((horizon+1)*2,1);-HH*lastline+KK];
% X=Abigx0+Bbig*U
% XredH*x<=XredK
% XredHbig*X<=XredKbig
% XredHbig*(Abigx0+Bbig*U)<=XredKbig
% XredHbig*Bbig*U<=XredKbig-XredHbig*Abigx0
  XredHbig=kron(eye(horizon+1),XredH);
  XredKbig=kron(ones(horizon+1,1),XredK);
  
  A=[A;XredHbig*Bbig;HH*Bbig(end-1:end,:)];
  b=[b;XredKbig-XredHbig*Abig*x0;-HH*lastline+KK];



    % f, the linear term incluedes also the refx and refu since we are
    % penalizing the distance from the reference. this latter can change
    % over the horizon
    % note that without preview the reference is the same over the entire
    % horizon while with preview it sees the future reference
    % in this last case, the state will start following the futur reference
    % in advance
    
    f=2*(x0'*Abig'*Qbig*Bbig);
    f=f';
    
    [U,FVAL,EXITFLAG]=quadprog(2*H,f,A,b);
    
    if EXITFLAG==-2||isempty(U)% in this cases the problem was not feasible
        flagtot=0;
        break
    else
        %from U we reconstruct the entire prediction of X
        X=Abig*x0+Bbig*U;
        % here below we reshape the prediction to have the states divided
        % on rows
        X_re=[reshape(X,2,horizon+1)'];
        % here we store in xx the sequence of closed loop states by
        % applying the first input only
        xx=[xx;X_re(2,:)];
        uu=[uu;U(1,:)];
    end
    
end


% plotting
t=0:Ts:simul_length*Ts-Ts;
subplot(3,1,1)
plot(t(1:size(xx,1)),xx(1:size(xx,1),1))
hold on
plot(t(1:size(xx,1)),repmat(xbounds,size(xx,1),1),'r--');
plot(t(1:size(xx,1)),repmat(-xbounds,size(xx,1),1),'r--');


grid on
title('X1 vs tiempo')
legend('x1 al aplicar MPC')
subplot(3,1,2)
plot(t(1:size(xx,1)),xx(1:size(xx,1),2))
hold on
plot(t(1:size(xx,1)),repmat(xbounds,size(xx,1),1),'r--');
plot(t(1:size(xx,1)),repmat(-xbounds,size(xx,1),1),'r--');

grid on
title('X2 vs tiempo')
legend('x2 al aplicar MPC')

subplot(3,1,3)
plot(t(1:size(uu,1)),uu(1:size(uu,1)));
hold on
plot(t(1:size(uu,1)),repmat(ubounds,size(uu,1),1),'r--');
plot(t(1:size(uu,1)),repmat(-ubounds,size(uu,1),1),'r--');

grid on
title('u vs tiempo')
legend('u al aplicar MPC')


