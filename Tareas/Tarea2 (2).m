clear all
close all
A=[4/3 -2/3;1 0];
B=[1 0]';
Q=[4/9+0.001 -2/3;-2/3 1.001];
R=0.001;

N=6;

Ap=zeros((N-1)*length(A),length(A));
Bp=zeros(length(Ap),(N-1));

Qp=zeros((N-1)*2,(N-1)*2);
Qp=blkdiag(Q,Q,Q,Q,Q);

Rp=zeros(N-1,N-1);
Rp=diag([R,R,R,R,R]);



%X=Axt|t+BU Vectorized notation

%Matrix A
for i=1:N-1
  %if i==1
    %Ap(i:i+1,:)=eye(length(A));
    
  %else
     
    Ap(i+(i-1):2*i,:)=A^(i-1);
    
  %end
  
end


%Bp(3:4,1)=(A^0)B;
%Bp(5:6,1)=A*B;
%Bp(7:8,1)=(A^2)*B;
%Bp(9:10,1)=(A^3)*B;

%Bp(5:6,2)=B;
%Bp(7:8,2)=A*B;
%Bp(9:10,2)=(A^2)*B;

%Bp(7:8,3)=B;
%Bp(9:10,3)=A*B;

%Bp(9:10,4)=B;

n=1;
%Matrix B
for i=1:N-2
    for j=n:N-2
       Bp(2*j+1:2*j+2,i)=A^(j-n)*B;
       
    end
    n=n+1;
end

