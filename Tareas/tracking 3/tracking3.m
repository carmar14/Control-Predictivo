clc
close all
clear


Ts=0.1; %sampling time
Ad=[1.7 -0.4;0.4 0.8]; % discrete time system xk+1=Axk+Buk
Bd=[0.2;0.05];          %                      yk=Cxk+Duk
Cd=eye(2);
Dd=0;

% initial condition
x0=[0.07501;-0.3327];%[-1,-1];%[0.07501;-0.3327];%[0.3;0.1];

%simulation length
simul_length=300;

% steady inputs
Us=[-0.8 -0.75 -0.65 -0.55 -0.45 -0.35];
% associated steady states
Xs=inv(eye(2)-Ad)*Bd*Us;


% constraints
ubounds=1; % bound on the absolute value of u
xbounds=2;

% MPC parameters
horizon=15;
flag_with_preview=1;
Q=diag([1;1]); % Q and R matrices of sum x'Qx+u'Ru
R=1; % please change here to see the effect of different values of Q and R

% LQR solution used for the terminal penalty and terminal constraint
[K,P]=dlqr(Ad,Bd,Q,R); % terminal penalty for stability

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

% Definition of the reference over time
refx=[];
refu=[];
for n=1:simul_length+horizon
    if n>=1 && n<=50
        refx=[refx Xs(:,1)];
        refu=[refu Us(1)];
    elseif n>=51 && n<=100
        refx=[refx Xs(: ,2)];
        refu=[refu Us(2)];             
    elseif n>=101 && n<=150
        refx=[refx Xs(:,3)];
        refu=[refu Us(3)];        
    elseif n>=151 && n<=200
        refx=[refx Xs(:,4)];
        refu=[refu Us(4)];        
    elseif n>=201 && n<=250
        refx=[refx Xs(:,5)];
        refu=[refu Us(5)];        
    elseif n>=251 && n<=simul_length+horizon
        refx=[refx Xs(:,6)];
        refu=[refu Us(6)];
    end    
end


xx(1,:)=x0';
uu=[];


for ii=1:simul_length-1
    
    % definition of the terminal set used at the end of each horizon
    if flag_with_preview==1
        [temp HH KK set_X]=fsetX(K,ubounds,xbounds,refx(:,ii+horizon),refu(:,ii+horizon));
    else
        [temp HH KK set_X]=fsetX(K,ubounds,xbounds,refx(:,ii),refu(:,ii));
    end
    
    x0=xx(ii,:)';
    
    % definition of the bounds on the input
    A=kron(eye(horizon),[1;-1]);
    b=ubounds*ones(horizon*2,1);
    
    % temp and lastline are used to impose the terminal constraint 
    temp=Abig*x0;
    lastline=temp(end-1:end,:);
    
    % the two lines below include the previous input constraints, the state
    % constraints expressed in terms of inputs and the terminal constraint
    A=[A;Bbig;-Bbig;HH*Bbig(end-1:end,:)];
    b=[b;-Abig*x0+xbounds*ones((horizon+1)*2,1);Abig*x0+xbounds*ones((horizon+1)*2,1);-HH*lastline+KK];
    
    % f, the linear term incluedes also the refx and refu since we are
    % penalizing the distance from the reference. this latter can change
    % over the horizon
    % note that without preview the reference is the same over the entire
    % horizon while with preview it sees the future reference
    % in this last case, the state will start following the futur reference
    % in advance
    
    if flag_with_preview==1
        f=2*(x0'*Abig'*Qbig*Bbig-reshape(refu(:,ii:ii+horizon-1),(horizon)*size(Bd,2),1)'*Rbig-reshape(refx(:,ii:ii+horizon),(horizon+1)*size(Ad,1),1)'*Qbig*Bbig);
    else
        f=2*(x0'*Abig'*Qbig*Bbig-reshape(repmat(refu(:,ii),horizon,1),(horizon)*size(Bd,2),1)'*Rbig-reshape(repmat(refx(:,ii),horizon+1,1),(horizon+1)*size(Ad,1),1)'*Qbig*Bbig);
    end
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
subplot(2,1,1)
plot(t(1:size(xx,1)),refx(1,1:size(xx,1)),'--',t(1:size(xx,1)),xx(1:size(xx,1),1))
grid on
if flag_with_preview==1
    title('X1 vs tiempo con preview')
else
    title('X1 vs tiempo sin preview')
end
legend('Referencia X1','x1 al aplicar MPC')
subplot(2,1,2)
plot(t(1:size(xx,1)),refx(2,1:size(xx,1)),'--',t(1:size(xx,1)),xx(1:size(xx,1),2))
grid on
if flag_with_preview==1
    title('X2 vs tiempo con preview')
else
    title('X2 vs tiempo sin preview')
end
legend('Referencia X2','x2 al aplicar MPC')




