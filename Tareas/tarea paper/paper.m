clear
close all
clc

A=[1 1; 0 1];
B=[0 0.5;1 0.5];
C=[1 0];
D=[0 0];

Q=eye(2);
R=eye(2);



% LQR solution used for the terminal penalty and terminal constraint
[K,P]=dlqr(A,B,Q,R); % terminal penalty for stability

T=100*P;

K=-K;
%Referencias
x0=[2;0];
xs=[4.95;0];
us=[0;0];
%Restricciones
xbounds=5;
ubounds=0.3;
%Encontrando el Kernel
matriz=[A-eye(size(A)) B];
Mtetha=null(matriz); %[xs1;xs2;u1s;u2s]
%Mtetha=[1 0 0 0;0 1 1 -2]';
Ntetha=C*Mtetha(1:2,:)+D*Mtetha(2:3,:);

%horizonte
N=3;
lambda=0.99;

%%%% Sistema aumentado
L=[-K eye(2)]*Mtetha;

Aw=[A+B*K B*L; zeros(2,2) eye(2)];
[temp HH KK set_w]=fsetW(K,lambda,L,ubounds,xbounds,Aw);

% set_w{end}.projection([1 2]).plot

% test_invariancia_tracking
% [temp HH KK set_X]=fsetX(K,ubounds,xbounds,xs,us)

% vectorial formulation
Qbig=kron(eye(N),Q);
Qbig=blkdiag(Qbig,P);
Q=Qbig;
Rbig=kron(eye(N),R);

Qbig=blkdiag(Qbig,Rbig);

Tbig=kron(eye(N+1),T);%repmat(T,size(Qbig)/2);

Mtethax=Mtetha(1:size(A),1:size(Mtetha,2));
Mtethau=Mtetha(size(A)+1:size(Mtetha,1),1:size(Mtetha,2));
% Mtethaxb=repmat(Mtethax',size(Q)/2);
% Mtethaub=repmat(Mtethau',size(Rbig)/2);
Mtethaxb=kron(eye(N+1),Mtethax');
Mtethaub=kron(eye(N),Mtethau');

MtxQ=-Mtethaxb'*Q;
MtuR=-Mtethaub'*Rbig;
last=Mtethaxb'*Q*Mtethaxb+[Mtethaub'*Rbig*Mtethaub zeros(6,2);zeros(2,8)]+Mtethaxb'*Tbig*Mtethaxb;
Qbbig=[Qbig [MtxQ;MtuR zeros(size(MtuR,1),2)];[MtxQ [MtuR zeros(size(MtuR,1),2)]'] last];
%termino lineal
f=[zeros(1,size(A,1)*(N+1)) zeros(1,size(B,2)*N) (-2*xs'*Mtethax'*T)*ones(2,size(A,1)*(N+1))];


%Matrices aumentadas

Abig=[];
Abig=kron(eye(N),A);
[m n]=size(A);
% Abig=[zeros(2,n+2); Abig];
A2=zeros((N+1)*m,(N+1)*m);
A2(m+1:end,1:end-n)=Abig;
Abig=A2;

Bbig=[];
Bbig=kron(eye(N),B);
[m n]=size(B);
% Abig=[zeros(2,n+2); Abig];
B=zeros((N+1)*m,(N+1)*n);
B(m+1:end,1:end-n)=Bbig;
Bbig=B;
Bbig=Bbig(:,1:end-n);

Hbig=[];
Hbig=eye(m);
H=zeros((N+1)*m,m);
H(1:m,1:m)=Hbig;
Hbig=H;
%cost function termino cuadratico
H=Qbbig;


%Restriciones desigualdad
matriz=zeros(size(A,1),((N+1)*size(A,1)));
matriz(:,end-size(A,1)+1:end)=eye(size(A,1)) ;

Ex=kron(eye(N+1),[1 0;-1 0;0 1;0 -1]);
Eu=kron(eye(N),[1 0;-1 0;0 1;0 -1]);
% Etethax=kron(eye(N+1),[1 0;-1 0;0 1;0 -1]);
% Etethau=kron(eye(N),[1 0;-1 0;0 1;0 -1]);
% Etetha=[Etethax];
% Ez=blkdiag([Ex;HH*matriz],Eu,Etethax);%Etethau);  %A para quadprog
Ez=blkdiag([Ex;HH*matriz],Eu,HH*matriz);
bu=ubounds*ones(N*4,1);
xu=xbounds*ones((N+1)*4,1);
% Fz=[xu;KK;bu;zeros(16,1)]; % b para quadprog duda para restriccion de desigualdad nueva
Fz=[xu;KK;bu;KK];
%Restriccions de igualdad
[m n]=size(Abig);
I=eye(size(Abig,2));
% Mtethab=kron(eye(2),Mtetha);
% Aeq=[[I-Abig; zeros(14,8)] [-Bbig; zeros(14,6)] [zeros(14,4);Mtethab] zeros(22,4)];%Mtethaxb [Mtethaub; zeros(2,6)]];
Mtethab=kron(eye(N+1),Mtetha);
Aeq=[[I-Abig; zeros(size(Mtethab))] [-Bbig; zeros(size(Mtethab,1),size(Bbig,2))] [zeros(size(Abig));Mtethab]];%Mtethaxb [Mtethaub; zeros(2,6)]];



% x0=[]

% [U,FVAL,EXITFLAG]=quadprog(2*H,f,Ez,Fz,Aeq,beq);

simul_length=50;
Tf=simul_length;
x1=[x0(1,1)];
x2=[x0(2,1)];
B=[0 0.5;1 0.5];
ss=[xs;us];
ssb=kron(ones(N+1,1),ss);
for ii=1:simul_length-1
    
    %     beq=[Hbig*x0; zeros(14,1)]; %duda para restriccion de igualdad nueva
%         beq=[Hbig*x0; zeros(16,1)];
    beq=[Hbig*x0; ssb];
    [U,FVAL,EXITFLAG]=quadprog(2*H,f,Ez,Fz,Aeq,beq);
    
    if EXITFLAG==-2||isempty(U)
        flagtot=0;
    else
        x_t1=A*x0+B*U(N*size(A,2)+3:N*size(A,2)+4); %Midiendo el nuevo estado a partir de la accion de control proveniente del problema de opt.
        x1=[x1 x_t1(1,1)];
        x2=[x2 x_t1(2,1)];
        x0=x_t1;
        
    end
end

% plotting
close all
figure
subplot(2,1,1)
Ts=1;
t=0:Ts:Tf*Ts-Ts;
hold on
plot(x1)%,'-',t,refx(1,1:end-horizon),'--')
grid on
% legend('Referencia X1 alt','x1 al aplicar MPC alt')
subplot(2,1,2)
plot(x2)%,'-k',t,refx(2,1:end-horizon),'-g')
grid on
% legend('Referencia X2 alt','x2 al aplicar MPC alt')

