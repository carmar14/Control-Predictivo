clear 
close all
clc

%Matrices A,B,Q
Ad=[4/3 -2/3; 1 0];
Bd=[1 0]';
Q=[4/9+0.001 -2/3; -2/3 1.001];
C=[-2/3 1];
R=0.001;
x0=[-0.5 1.3]';
ubounds=1; % Valores absolutos de las restricciones de u
xbounds=2;

N=30;

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
%Rx*X<=k, restricciones para los estados
%Ru*U<=Ku, restriciones para la entrada

A=kron(eye(N),[1;-1]);
b=ubounds*ones(N*2,1);

A=[A;Bbig;-Bbig];
b=[b;-Abig*x0+xbounds*ones(N*2,1);Abig*x0+xbounds*ones(N*2,1)];
f=Bbig'*Qbig*Abig*x0;

[U,FVAL,EXITFLAG]=quadprog(2*H,2*f,A,b);
Ts=0.1;
if EXITFLAG==-2||isempty(U)
    flagtot=0;
else
    X=Abig*x0+Bbig*U;
    X_re=[reshape(X,2,N)'];
    
    figure(1)
    
    hold on
    subplot(1,2,1);plot(0:Ts:Ts*(N-1),X_re(:,1),'b');hold on;plot(0:Ts:Ts*(N-1),X_re(:,1),'bo');plot(0:Ts:Ts*(N-1),X_re(:,2),'r');plot(0:Ts:Ts*(N-1),X_re(:,2),'ro');grid on
    subplot(1,2,2);plot(0:Ts:Ts*(N-1),U,'k');hold on;plot(0:Ts:Ts*(N-1),U,'ko');grid on
    
    
end

% Control 
Tf=50;
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
        x_t1=Ad*x0+Bd*U(1);
        x1=[x1 x_t1(1,1)];
        x2=[x2 x_t1(2,1)];
        x0=x_t1;
    end
    
    clc
%     X_re=[reshape(X,2,N)'];
end

figure
t=0:Ts:Tf*Ts;%-Ts;
hold on
plot(t,x1,'-o')
plot(t,x2,'-ko')  
grid on
legend('x1 state','x2 state')
    

%Forma alterna



% N=6;   %5+1 por la programacion de matlab
% 
% %Nuevas matrices
% An=zeros(length(A)*(N-1),2);
% j=0;
% for i=1:(length(An))/2
% %     if i==1
% %         An(i:2,:)=eye(2);
% %     else
%         An(i+j:i+j+1,:)=A^(i-1);
%         j=j+1;
% %     end
% end
% 
% Bn=zeros(length(A)*(N-2),5);
% 
% % j=0
% % for i=1:N-2
% % %     if i==1
% % %         Bn(i,:)=0;
% % %     else
% % %         
% % %     for j=1:N-2
% %     
% %         Bn(i+j+1:i+j+2,1)=A^(i-1)*B
% %         A^(i-1)*B
% %         j=j+1;
% % %     end
% %     
% % end
% % % Bn
% n=1;
% for i=1:N-2
%     for j=n:N-2
%        Bn(2*j+1:2*j+2,i)=A^(j-n)*B;
%        
%     end
%     n=n+1;
% end
% 
% Qn=zeros((N-1)*2,(N-1)*2);
% Qn=blkdiag(Q,Q,Q,Q,Q);
% 
% Rn=zeros(N-1,N-1);
% Rn=diag([R,R,R,R,R]);
% Flu=-1;
% Fuu=1;
% 
% H=2*Bn'*Qn*Bn+R;
% f=(Bn'*Qn*An)'; %xt|t????
% 
% 
