clc
close all
clear

Ts=0.1; %sampling time
Ad=[1.7 -0.4;0.4 0.8]; % discrete time system xk+1=Axk+Buk
Bd=[0.2;0.05];          %                      yk=Cxk+Duk
Cd=eye(2);
Dd=0;

horizon=15;
Us=[-0.8 -0.5 0 -0.4 0.2 -0.3 0.2 0.5 -0.2 0.4];

Xs=inv(eye(2)-Ad)*Bd*Us;

ubounds=1; % bound on the absolute value of u
xbounds=2;

Q=diag([1;1]); % Q and R matrices of sum x'Qx+u'Ru
R=1; % please change here to see the effect of different values of Q and R
[K,S,E] = dlqr(Ad,Bd,Q,R);

[temp HH KK set_X]=fsetX(K,ubounds,xbounds,Xs(:,1),Us(1));



Qbig=kron(eye(horizon),Q);

[~,P]=dlqr(Ad,Bd,Q,R); % terminal penalty for stability
% the terminal constraint is not implemented yet
% P=S;
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

% constraints on the inputs - no state constraints for the moment. Please
% add them


x0=[0.155 ;0.602];%[0.3;0.1];[0.07501 ;-0.3327]
simul_length=100;
x2=[];
un2=[];
x1d=[];
x2d=[];
un=[];
refx=[];
refu=[];

for i=1:simul_length
    if i>=1 && i<=10
        refx=[refx Xs(:,1)];
        refu=[refu Us(1)];
    elseif i>10 && i<=20
        refx=[refx Xs(:,2)];
        refu=[refu Us(2)];
    elseif i>20 && i<=30
        refx=[refx Xs(:,3)];
        refu=[refu Us(3)];
    elseif i>30 && i<=40
        refx=[refx Xs(:,4)];
        refu=[refu Us(4)];
    elseif i>40 && i<=50
        refx=[refx Xs(:,5)];
        refu=[refu Us(5)];
    elseif i>50 && i<=60
        refx=[refx Xs(:,6)];
        refu=[refu Us(6)];
    elseif i>60 && i<=70
        refx=[refx Xs(:,7)];
        refu=[refu Us(7)];
    elseif i>70 && i<=80
        refx=[refx Xs(:,8)];
        refu=[refu Us(8)];
    elseif i>80 && i<=90
        refx=[refx Xs(:,9)];
        refu=[refu Us(9)];
    else
        refx=[refx Xs(:,10)];
        refu=[refu Us(10)];
    end

end

xsbig=repmat(refx,horizon+1,1);
usbig=repmat(refu,horizon,1);
xx(1,:)=x0';
uu=[];
o=1;

for ii=1:simul_length

    if mod(ii,10)==0
    [temp HH KK set_X]=fsetX(K,ubounds,xbounds,Xs(:,o),Us(o));
    o=o+1;
    end
    x0=xx(ii,:)';

    A=kron(eye(horizon),[1;-1]);
    b=ubounds*ones(horizon*2,1);

    temp=Abig*x0;
    lastline=temp(end-1:end,:);

    A=[A;Bbig;-Bbig;HH*Bbig(end-1:end,:)];
    b=[b;-Abig*x0+xbounds*ones((horizon+1)*2,1);Abig*x0+xbounds*ones((horizon+1)*2,1);-HH*lastline+KK];

    if ii==simul_length
        break
    end
    f=2*(x0'*Abig'*Qbig*Bbig-usbig(:,ii+1)'*Rbig-xsbig(:,ii+1)'*Qbig*Bbig);%((x0'*Abig'-xs')*Qbig*Bbig-us'*Rbig);%Bbig'*Qbig*Abig*x0;
    f=f'; 

    [U,FVAL,EXITFLAG]=quadprog(2*H,f,A,b);

    if EXITFLAG==-2||isempty(U)
        flagtot=0;
        break
    else

    X=Abig*x0+Bbig*U;
    
    X_re=[x0';reshape(X,2,horizon+1)'];
    xx=[xx;X_re(end,:)];
    uu=[uu;U(1,:)];
    end
end

% close all
figure

t=0:Ts:simul_length*Ts-Ts;
subplot(2,1,1)
plot(t,refx(1,:),'--',t,xx(1:length(t),1))
grid on
title('X1 vs tiempo')
legend('Referencia X1','x1 al aplicar MPC')
subplot(2,1,2)
plot(t,refx(2,:),'--',t,xx(1:length(t),2))
grid on
title('X2 vs tiempo')
legend('Referencia X2','x2 al aplicar MPC')

% figure
% plot(t,refu,'--',t(1:length(uu)),uu)
% title('u vs tiempo')
% legend('Referencia u','u del MPC')

