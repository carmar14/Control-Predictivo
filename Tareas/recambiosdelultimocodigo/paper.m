clear
close all
clc

A=[1 1; 0 1];
B=[0 0.5;1 0.5];
C=[1 0];
D=[0 0];

Qi=eye(2);
R=eye(2);



% LQR solution used for the terminal penalty and terminal constraint
[K,P]=dlqr(A,B,Qi,R); % terminal penalty for stability

T=100*P;

K=-K;
%Referencias
x0=[2;0];
x0=[0;0];

%Restricciones
xbounds=5;
ubounds=0.3;
%Encontrando el Kernel
%matriz=[A-eye(size(A)) B];
%Mtetha=null(matriz); %[xs1;xs2;u1s;u2s]
%Mtetha=[1 0 0 0;0 1 1 -2]';

%%% Comentarios: no. Ntheta no se calcula asi'. Si uno necesita tambien
%%% Ntetha hay que resolver (3) asi'
Matrix_theta=null([A-eye(2) B zeros(size(A,2),1);C D -eye(size(C,1))]);
Mtetha=Matrix_theta(1:size(A,2)+size(B,2),:);

%%%Ntheta=Matrix_theta(size(A,2)+size(B,2)+1:end,:);
%Ntetha=C*Mtetha(1:2,:)+D*Mtetha(2:3,:);

%horizonte
N=3;
lambda=0.99;

%%%% Sistema aumentado
L=[-K eye(size(B,2))]*Mtetha;

Aw=[A+B*K B*L; zeros(2,2) eye(2)];
% [temp HH KK temp2 HH2 KK2 set_w]=fsetW(K,lambda,L,ubounds,xbounds,Aw);
[temp HH KK temp2 HH2 KK2 set_w]=fsetW(K,lambda,L,ubounds,xbounds,Aw,Mtetha);

theta=[-1;3];
xs=Mtetha(1:2,:)*theta;
us=Mtetha(3:4,:)*theta;

xs=[4.95; 0];
% set_w{end}.projection([1 2]).plot
figure
plot(set_w{end}.projection([1 2]))

% test_invariancia_tracking
% [temp HH KK set_X]=fsetX(K,ubounds,xbounds,xs,us)

% vectorial formulation
Qbig=kron(eye(N),Qi);
Qbig=blkdiag(Qbig,P);
Q=Qbig;
Rbig=kron(eye(N),R);

Qbig=blkdiag(Qbig,Rbig);

Tbig=kron(eye(N+1),T);%repmat(T,size(Qbig)/2);

Mtethax=Mtetha(1:size(A),1:size(Mtetha,2));
Mtethau=Mtetha(size(A)+1:size(Mtetha,1),1:size(Mtetha,2));
% Mtethaxb=repmat(Mtethax',size(Q)/2);
% Mtethaub=repmat(Mtethau',size(Rbig)/2);
% Mtethaxb=kron(eye(N+1),Mtethax');
% Mtethaub=kron(eye(N),Mtethau');
Mtethaxb=Mtethax'*Qi;
Mtethaub=Mtethau'*R;

MtxQ=repmat(-Mtethaxb,1,N);
MtuR=repmat(-Mtethaub,1,N);
MtxQ2=repmat(-Mtethaxb,N+1,1);
MtuR2=repmat(-Mtethaub,N,1);
%last=(Mtetha(1:2,:)'*Qi)*Mtetha(1:2,:)*N+(Mtetha(3:4,:)'*R)*Mtetha(3:4,:)*N+(Mtetha(1:2,:)'*P)*Mtetha(1:2,:)+(Mtetha(1:2,:)'*T)*Mtetha(1:2,:);%Mtethax'*Qi*Mtethax+Mtethau'*R*Mtethau+Mtethax'*T*Mtethax;
last=(Mtethaxb*Mtethax)*N+(Mtethaub*Mtethau)*N+Mtethax'*P*Mtethax+Mtethax'*T*Mtethax;
% Qbbig=[Qbig [MtxQ2; MtuR2];MtxQ MtuR last];
Qbbig=[Qbig [MtxQ -Mtethax'*P MtuR]';MtxQ -Mtethax'*P MtuR last];
%termino lineal
%f=[zeros(1,size(A,1)*(N+1)) zeros(1,size(B,2)*N) (-2*xs'*Mtethax'*T)*ones(2,size(A,1)*(N+1))];
f=[zeros(1,size(A,1)*(N+1)) zeros(1,size(B,2)*N) -2*xs'*T*Mtethax ];%-2*xs'*Mtethax'*T ];
%%% Comentarios
%%% He mirado las dimensiones de Qbbig y f y hay algo que no entiendo
%%% si N=3, y tenemos 2 estados, 2 entradas y theta de dimension 2,
%%% entonces tenemos un total de (N+1)*2+N*2+2=16 mientras aqui' hay 22..


% keyboard

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
B2=zeros((N+1)*m,(N+1)*n);
B2(m+1:end,1:end-n)=Bbig;
Bbig=B2;
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
HHx=HH(:,1:2);
HHtetha=HH(:,3:end);
% Ez=blkdiag([Ex;HHx*matriz],Eu,HHtetha*eye(2));

Ez=blkdiag(Ex,Eu);
Ez=[Ez zeros(size(Ez,1),size(Mtetha,2))];

Ez=[Ez;[zeros(size(HHx,1),size(A,1)*N) HHx zeros(size(HHx,1),size(B,2)*N) HHtetha]];

%%% Comentarios: Cuidado se necesitan las restricciones sobre X y Theta y no solo la
%%% proyeccion en X!


bu=ubounds*ones(N*4,1);
xu=xbounds*ones((N+1)*4,1);
% Fz=[xu;KK;bu;zeros(16,1)]; % b para quadprog duda para restriccion de desigualdad nueva
% KKx=KK(:,1:2);
% KKtetha=KK(:,3:end);
Fz=[xu;bu;KK];
%Restriccions de igualdad
[m n]=size(Abig);
I=eye(size(Abig,2));
% Mtethab=kron(eye(2),Mtetha);
% Aeq=[[I-Abig; zeros(14,8)] [-Bbig; zeros(14,6)] [zeros(14,4);Mtethab] zeros(22,4)];%Mtethaxb [Mtethaub; zeros(2,6)]];
% Mtethab=kron(eye(N+1),Mtetha);
% Aeq=[[I-Abig; zeros(size(Mtethab))] [-Bbig; zeros(size(Mtethab,1),size(Bbig,2))] [zeros(size(Abig));Mtethab]];%Mtethaxb [Mtethaub; zeros(2,6)]];
Mtethab=[Mtetha ; Mtetha];
%Aeq=[I-Abig -Bbig Mtethab];
Aeq=[I-Abig -Bbig zeros(size(Abig,1),size(Mtetha,2))];

% x0=[]

% [U,FVAL,EXITFLAG]=quadprog(2*H,f,Ez,Fz,Aeq,beq);

simul_length=90;
Tf=simul_length;
x1=[x0(1,1)];
x2=[x0(2,1)];
%B=[0 0.5;1 0.5];
% ss=[xs;us];
% ssb=kron(ones(N+1,1),ss);
k=1;
%%% Comentarios: no entiend bien porque en beq entra ssb en la forma que
%%% pusiste aqui. En principio entra multiplicando -2*Theta'*Mthetax'*T*xs.
%%% Asi que no hace falta us y no deberia aparecer como otros elementos de
%%% beq
xs=[];
for i=1:Tf
    if (i>=1 && i<=30)
        xs(:,i)=[4.95;0];
    end
    if (i>30 && i<=60)
        xs(:,i)=[-5.5;0];
    end
    if(i>60)
        xs(:,i)=[2;0];
    end
    
end
%%Estados artificiales
xss=[];



for ii=1:simul_length-1
    
    %     beq=[Hbig*x0; zeros(14,1)]; %duda para restriccion de igualdad nueva
    %         beq=[Hbig*x0; zeros(16,1)];
    beq=[Hbig*x0] ; %ones((N+1)*size(A,1),1)];
    f=[zeros(1,size(A,1)*(N+1)) zeros(1,size(B,2)*N) -2*xs(:,ii)'*T*Mtethax ];

    [U,FVAL,EXITFLAG]=quadprog(2*H,f,Ez,Fz,Aeq,beq);
    
    if EXITFLAG==-2||isempty(U)
        flagtot=0;
    else
        
        x_t1=A*x0+B*U(N*size(A,2)+3:N*size(A,2)+4); %Midiendo el nuevo estado a partir de la accion de control proveniente del problema de opt.
        x1=[x1 x_t1(1,1)];
        x2=[x2 x_t1(2,1)];
        x0=x_t1;
        uc(:,k)=U(N*size(A,2)+3:N*size(A,2)+4);
        k=k+1;
        xss(:,k)=Mtetha*U(end-1:end,:);
    end
end

% plotting
%close all
figure
subplot(2,1,1)
Ts=1;
t=0:Ts:Tf*Ts-Ts;

plot(x1)%,'-',t,refx(1,1:end-horizon),'--')
hold on
plot(xss(1,:),'--');
plot(xs(1,:),'--k');
legend('Estado1 Real','Estado1 Artificial','Referencia');
ylim([-6 6]);
grid on
% legend('Referencia X1 alt','x1 al aplicar MPC alt')
subplot(2,1,2)
plot(x2)%,'-k',t,refx(2,1:end-horizon),'-g')
grid on
hold on
plot(xss(2,:),'--');
plot(xs(2,:),'--k');
legend('Estado2 Real','Estado2 Artificial','Referencia');
% legend('Referencia X2 alt','x2 al aplicar MPC alt')

