clc
close all
clear all

Ts=0.1; %sampling time
Ad=[1.7 -0.4;0.4 0.8]; % discrete time system xk+1=Axk+Buk
Bd=[0.2;0.05];          %                      yk=Cxk+Duk
Cd=eye(2);
Dd=0;

x0=[0.3;-0.1]; % initial condition of the system

ubounds=1; % bound on the absolute value of u
xbounds=2;

Q=diag([1;1]); % Q and R matrices of sum x'Qx+u'Ru
R=1; % please change here to see the effect of different values of Q and R
[K,S,E] = dlqr(Ad,Bd,Q,R);

% -ubounds<=-Kx<=ubounds -xbounds<=x<=xbounds
setXU=Polyhedron('A',[-K;K;eye(2);-eye(2)],'b',[ubounds;ubounds;xbounds;xbounds;xbounds;xbounds]);
set_X{1}=setXU;

hold on
for i=1:100
    if i>1 && abs(volume(set_X{i-1})-volume(set_X{i}))<1e-3
        break
    end
    volume(set_X{i})
    try
        set_X{i+1}=inv(Ad-Bd*K)*set_X{i}&set_X{i};
            plot(set_X{i+1},'Color','r');
%             keyboard
    catch
        break
    end
end

temp=set_X{end}.H;
HH=temp(:,1:end-1);
KK=temp(:,end);

horizon=input('Which horizon? '); % the user can choose the horizon


%%
% % sys = ss(Ad,Bd,Cd,Dd,Ts);
% % model = LTISystem(sys);
% % model.x.min = -xbounds*ones(2,1);
% % model.x.max = xbounds*ones(2,1);
% % model.u.min = -ubounds;
% % model.u.max = ubounds;
% % model.x.penalty = QuadFunction(Q);
% % model.u.penalty = QuadFunction(R);
% % Tset = model.LQRSet;
% % PN = model.LQRPenalty;
% % 
% % model.x.with('terminalSet');
% % model.x.terminalSet = Tset;
% % model.x.with('terminalPenalty');
% % model.x.terminalPenalty = PN;
% % mpc = MPCController(model, horizon)
% % expmpc = mpc.toExplicit();

%%

figure(2)
% plot(setXU,'Color','r');
hold on
plot(set_X{end},'Color','y','FaceColor','none');
% plot(Tset,'Color','c');

% below we construct all the matrices according to what described in the
% scripts

Qbig=kron(eye(horizon-1),Q);

[~,P]=dlqr(Ad,Bd,Q,R); % terminal penalty for stability
% the terminal constraint is not implemented yet

Qbig=blkdiag(Qbig,P);
Rbig=kron(eye(horizon),R);
Abig=[];
Bbig=0*kron(eye(horizon),Bd);

for i=1:horizon
    Abig=[Abig;Ad^i];
    Bbig=Bbig+kron(diag(ones(horizon-i+1,1),1-i),Ad^(i-1)*Bd);
end



%cost function - we assume the presence of the quadratic term only

H=Bbig'*Qbig*Bbig+Rbig;

% constraints on the inputs - no state constraints for the moment. Please
% add them


for jj=1:1000
    clear xx uu 
    x0=2*xbounds*rand(2,1)-xbounds;
    
    xx(1,:)=x0';
    uu=[];
    simul_length=30; % length of the simulation. Here we apply the receding horizon approach
    % at each time step we solve the optimization over the
    % horizon. Then we apply the first input only and then we
    % re-optimize
    
    flagtot=1;
    
    for ii=1:simul_length
        x0=xx(ii,:)';
        
        A=kron(eye(horizon),[1;-1]);
        b=ubounds*ones(horizon*2,1);
        
        temp=Abig*x0;
        lastline=temp(end-1:end,:);
        
        
        
        A=[A;Bbig;-Bbig;HH*Bbig(end-1:end,:)];
        b=[b;-Abig*x0+xbounds*ones(horizon*2,1);Abig*x0+xbounds*ones(horizon*2,1);-HH*lastline+KK];
        
        
        f=Bbig'*Qbig*Abig*x0;
        
        
        [U,FVAL,EXITFLAG]=quadprog(2*H,2*f,A,b);
        
        if EXITFLAG==-2||isempty(U)
            flagtot=0;
            break
        else
            
        X=Abig*x0+Bbig*U;
        
        X_re=[x0';reshape(X,2,horizon)'];
        xx=[xx;X_re(end,:)];
        uu=[uu;U(1,:)];
            ;
        end
        
        
    end
    if not(flagtot==0)
    figure(1)
    
    hold on
    subplot(1,2,1);plot(0:Ts:Ts*(simul_length),xx(:,1),'b');hold on;plot(0:Ts:Ts*(simul_length),xx(:,1),'bo');plot(0:Ts:Ts*(simul_length),xx(:,2),'r');plot(0:Ts:Ts*(simul_length),xx(:,2),'ro');grid on
    subplot(1,2,2);plot(0:Ts:Ts*(simul_length-1),uu,'k');hold on;plot(0:Ts:Ts*(simul_length-1),uu,'ko');grid on
    end
    
    figure(2)
    if flagtot==1
        col='b';
        plot(xx(:,1),xx(:,2),col)
        hold on
        plot(xx(:,1),xx(:,2),strcat(col,'o'),'LineWidth',2)
        
    else
        col='r';
        plot(xx(1,1),xx(1,2),col)
        hold on
        plot(xx(1,1),xx(1,2),strcat(col,'o'),'LineWidth',2)
        
    end
    
    
    axis([-xbounds xbounds -xbounds xbounds])
    box on
    try
    close(1)
    catch
        ;
    end
end

figure(2)

plot(expmpc.partition)