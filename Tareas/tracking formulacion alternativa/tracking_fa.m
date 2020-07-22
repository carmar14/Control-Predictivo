clc
close all
clear


Ts=0.1; %sampling time
Ad=[1.7 -0.4;0.4 0.8]; % discrete time system xk+1=Axk+Buk
Bd=[0.2;0.05];          %                      yk=Cxk+Duk
Cd=eye(2);
Dd=0;
Q=diag([1;1]); % Q and R matrices of sum x'Qx+u'Ru
R=1; 

% MPC parameters
horizon=25;
flag_with_preview=0;


% initial condition
% x0=[0.07501;-0.3327]
x0=[0.07501;-0.3327]; 
% x0=[0.65;-0.33];%[-1,-1];%[0.07501;-0.3327];%[0.3;0.1];

%simulation length
simul_length=300;

% steady inputs
% us=0.8;

Us=[-0.6 -0.53 -0.27 0.2 0.38 0.52];
%Us=[0 0 0 0 0 0];
% associated steady states
%%Xs=inv(eye(2)-Ad)*Bd*Us;
Xs=inv(eye(2)-Ad)*Bd*Us;
% xsbig=repmat(xs,horizon+1,1);
% Xs=inv(eye(2)-Ad)*Bd*Us;

% constraints
ubounds=1; % bound on the absolute value of u
xbounds=2;


% Q=diag([1;1]); % Q and R matrices of sum x'Qx+u'Ru
% R=1; % please change here to see the effect of different values of Q and R

% LQR solution used for the terminal penalty and terminal constraint
[K,P]=dlqr(Ad,Bd,Q,R); % terminal penalty for stability

%[temp HH KK set_X]=fsetX(K,ubounds,xbounds,xs,us);

%[~,P]=dlqr(Ad,Bd,Q,R);
% vectorial formulation
Qbig=kron(eye(horizon),Q);
Qbig=blkdiag(Qbig,P);
Rbig=kron(eye(horizon),R);


Qbig=blkdiag(Qbig,Rbig);

Abig=[];
Abig=kron(eye(horizon),Ad);
[m n]=size(Ad);
% Abig=[zeros(2,n+2); Abig];
A=zeros((horizon+1)*m,(horizon+1)*m); 
A(m+1:end,1:end-n)=Abig;
Abig=A;

Bbig=[];
Bbig=kron(eye(horizon),Bd);
[m n]=size(Bd);
% Abig=[zeros(2,n+2); Abig];
B=zeros((horizon+1)*m,(horizon+1)*n); 
B(m+1:end,1:end-n)=Bbig;
Bbig=B;
Bbig=Bbig(:,1:end-n);

Hbig=[];
Hbig=eye(m);
H=zeros((horizon+1)*m,m);
H(1:m,1:m)=Hbig;
Hbig=[H ]; %zeros(size(Bbig,2),size(Ad,2))];

%cost function - we assume the presence of the quadratic term only

H=Qbig;


x1=[x0(1,1)];
x2=[x0(2,1)];
xx(1,:)=x0';
uu=[];
Tf=simul_length;

%f=zeros(3,horizon);

% Definition of the reference over time
refx=[];
refu=[];
refz=[];

for n=1:simul_length+horizon
    if n>=1 && n<=50
        refx=[refx Xs(:,1)];
        refu=[refu Us(1)];
        refz=[refx;refu];
    elseif n>=51 && n<=100
        refx=[refx Xs(: ,2)];
        refu=[refu Us(2)];
        refz=[refx;refu];
    elseif n>=101 && n<=150
        refx=[refx Xs(:,3)];
        refu=[refu Us(3)];
        refz=[refx;refu];
    elseif n>=151 && n<=200
        refx=[refx Xs(:,4)];
        refu=[refu Us(4)];
        refz=[refx;refu];
    elseif n>=201 && n<=250
        refx=[refx Xs(:,5)];
        refu=[refu Us(5)];
        refz=[refx;refu];
    elseif n>=251 && n<=simul_length+horizon
        refx=[refx Xs(:,6)];
        refu=[refu Us(6)];
        refz=[refx;refu];
    end    
end



for ii=1:simul_length-1
    
    % definition of the terminal set used at the end of each horizon
    if flag_with_preview==1
        [temp2 HH KK set_X]=fsetX(K,ubounds,xbounds,refx(:,ii+horizon),refu(:,ii+horizon));
    else
        [temp2 HH KK set_X]=fsetX(K,ubounds,xbounds,refx(:,ii),refu(:,ii));
    end
    

%      x0big=[];
%      for n=1:horizon+1
%          x0big=[x0big; x0];
%      end
%      
%     
%     % temp and lastline are used to impose the terminal constraint 
%     temp=Abig*x0big;
%     lastline=temp(end-1:end,:);
     
    %%Restricciones desigualdad
    %ExB=[Bbig;-Bbig];
    A=kron(eye(horizon),[1;-1]);
    Eu=A;
    
    
    matriz=zeros(size(Ad,1),((horizon+1)*size(Ad,1)));
    matriz(:,end-size(Ad,1)+1:end)=eye(size(Ad,1)) ;
    
    %% Faltaba algo aqui
    Ex=kron(eye(horizon+1),[1 0;-1 0;0 1;0 -1]);
    % Ex=[Ex;Ex];
    Ez=blkdiag([Ex;HH*matriz],Eu);  %A para quadprog
    %% error de dimension
    bu=ubounds*ones(horizon*2,1);
    xu=xbounds*ones((horizon+1)*4,1);
    %% no es un blockdiag aqui!
    Fz=[xu;KK;bu]; % b para quadprog

    %Restriccions de igualdad
    [m n]=size(Abig);
    I=eye(size(Abig,2));
    Aeq=[I-Abig -Bbig];
    beq=Hbig;
    
    if flag_with_preview==1
        f=-2*[reshape(refx(:,ii:ii+horizon),(horizon+1)*size(Ad,1),1);...
        reshape(refu(:,ii:ii+horizon-1),(horizon)*size(Bd,2),1)]'*Qbig;
    else
        f=-2*[reshape(repmat(refx(:,ii),horizon+1,1),(horizon+1)*size(Ad,1),1);...
        reshape(repmat(refu(:,ii),horizon,1),(horizon)*size(Bd,2),1)]'*Qbig;
    end
    
    
    
    f=f';
   
    [U,FVAL,EXITFLAG]=quadprog(2*H,f,Ez,Fz,Aeq,beq*x0);

    if EXITFLAG==-2||isempty(U)
        flagtot=0;
    else
        x_t1=Ad*x0+Bd*U(horizon*size(Ad,2)+3); %Midiendo el nuevo estado a partir de la accion de control proveniente del problema de opt.
        x1=[x1 x_t1(1,1)];
        x2=[x2 x_t1(2,1)];
        x0=x_t1;

    end
end


% plotting
close all
figure
subplot(2,1,1)
t=0:Ts:Tf*Ts-Ts;
hold on
plot(t,x1,'-',t,refx(1,1:end-horizon),'--')
grid on
legend('Referencia X1 alt','x1 al aplicar MPC alt')
subplot(2,1,2)
plot(t,x2,'-k',t,refx(2,1:end-horizon),'-g')  
grid on
legend('Referencia X2 alt','x2 al aplicar MPC alt')






