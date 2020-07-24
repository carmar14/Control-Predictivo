clear all
clc

A=[1 0.1 0;0.1 1 0.1;0 0.1 1];
Ac=A-diag(diag(A));
B=eye(3);



F{1}=[1;-1];% <=[1;1]
F{2}=[1;-1];% <=[1;1]
F{3}=[1;-1];% <=[1;1]
Xi{1}=1;
Xi{2}=1;
Xi{3}=1;
H{1}=[1;-1];% <=[1;1]
H{2}=[1;-1];% <=[1;1]
H{3}=[1;-1];% <=[1;1]


for i=1:3
    
    As{i}=A(i,i);
    N{i}=find(not(Ac(i,:)==0));
    Bs{i}=B(i,i);
    subsys{i}.mualfa=1;
    subsys{i}.mubeta=1;
    subsys{i}.Ad=As{i};
    subsys{i}.B=Bs{i};
    subsys{i}.Aij=A(i,N{i});
    subsys{i}.Fs=F{i};
    subsys{i}.Fc=F(N{i});
    subsys{i}.Xis=Xi{i};
    subsys{i}.Xic=Xi(N{i});
    subsys{i}.Xishat=Xi{i};
    subsys{i}.Xichat=Xi(N{i});
    subsys{i}.H=H{i};
    subsys{i}.w=Polyhedron('A',[],'b',[]);
    for j=1:length(N{i})
        subsys{i}.w=subsys{i}.w+subsys{i}.Aij(1,j)*Polyhedron('A',subsys{i}.Fc{j},'b',ones(size(subsys{i}.Fc{j},1),1));
    end
end
for i=1:3
    if isempty(null(ctrb(As{i},Bs{i})))
        disp('OK')
    else
        disp('Not OK');
        keyboard
    end
end
for i=1:3
    [var{i},fval{i},exitflag{i}] = fmincon(@(var) myfun(var,subsys{i}),[1;1;0.01],[],[],[],[],[1e-3;1e-3;eps],[1e3;1e3;1e-1],@(var) mycon(var,subsys{i}));
    [c{i}] = mycon(var{i},subsys{i});
end



function [K,P]=getKP(Qvec,Rvec,Asub,Bsub)
[K,P] = dlqr(Asub,Bsub,Qvec,Rvec);
K=-K;
end

function f=myfun(var,subsyst)%Asub y Bsub son parametros


Qvec=var(1:size(subsyst.Ad,1));
Rvec=var(size(subsyst.Ad,1)+1:size(subsyst.Ad,1)+size(subsyst.B,2));
delta=var(end);
[K,P]=getKP(Qvec,Rvec,subsyst.Ad,subsyst.B);
F=subsyst.Ad+subsyst.B*K;
Z=epsilon_mRPI(F,subsyst.w,delta);
[HH,KK]=Z.doubleHK;%Hx<=K
mRPIset=Polyhedron('A',HH,'b',KK);


alfa=0;
for j=1:size(subsyst.Aij,2)
    k=0;
%     disp(sprintf('Valor de j %d',j));
    while 1
        temp=norm(subsyst.Fs*F^k*subsyst.Aij(1,j)*pinv(subsyst.Fc{j}),Inf);
%         temp
        alfa=alfa+temp;
        if temp<1e-5
            break
        end        
        k=k+1;
    end
end





%%%%%%%%%
for r=1:size(subsyst.H,1)
    [variable, coste]=linprog(-subsyst.H(r)'*K,HH,KK);
    lhatv(r)=-coste;
end
beta=max(lhatv);

for r=1:size(subsyst.Fs,1)
    Lhatir(r)=1;
    for j=1:size(subsyst.Aij,2)
        k=0;
%         disp(sprintf('Valor de j %d',j));
        while 1
            temp=norm(subsyst.Fs(r,:)'*F^k*subsyst.Aij(1,j)*(subsyst.Xic{j}),Inf);
%             temp
            Lhatir(r)=Lhatir(r)-temp;
            if temp<1e-5
                break
            end
            k=k+1;
        end
    end
    Lhatir(r)=Lhatir(r)/norm(subsyst.Fs(r,:)'*subsyst.Xishat,Inf);
    diff_L(r)=Lhatir(r)-(norm(subsyst.Fs(r,:)',Inf)*delta/norm(subsyst.Fs(r,:)'*subsyst.Xishat,Inf));
end
Lhat=min(diff_L);
%%%%%%%%%
% lhatv
% Lhat
% mRPIset
f=subsyst.mualfa*alfa+subsyst.mubeta*beta;

end

function [c,ceq] = mycon(var,subsyst)


Qvec=var(1:size(subsyst.Ad,1));
Rvec=var(size(subsyst.Ad,1)+1:size(subsyst.Ad,1)+size(subsyst.B,2));
delta=var(end);
[K,P]=getKP(Qvec,Rvec,subsyst.Ad,subsyst.B);
F=subsyst.Ad+subsyst.B*K;
Z=epsilon_mRPI(F,subsyst.w,delta);
[HH,KK]=Z.doubleHK;%Hx<=K
mRPIset=Polyhedron('A',HH,'b',KK);


alfa=0;
for j=1:size(subsyst.Aij,2)
    k=0;
%     disp(sprintf('Valor de j %d',j));
    while 1
        temp=norm(subsyst.Fs*F^k*subsyst.Aij(1,j)*pinv(subsyst.Fc{j}),Inf);
%         temp
        alfa=alfa+temp;
        if temp<1e-5
            break
        end        
        k=k+1;
    end
end





%%%%%%%%%
for r=1:size(subsyst.H,1)
    [variable, coste]=linprog(-subsyst.H(r)'*K,HH,KK);
    lhatv(r)=-coste;
end
beta=max(lhatv);

for r=1:size(subsyst.Fs,1)
    Lhatir(r)=1;
    for j=1:size(subsyst.Aij,2)
        k=0;
%         disp(sprintf('Valor de j %d',j));
        while 1
            temp=norm(subsyst.Fs(r,:)'*F^k*subsyst.Aij(1,j)*(subsyst.Xic{j}),Inf);
%             temp
            Lhatir(r)=Lhatir(r)-temp;
            if temp<1e-5
                break
            end
            k=k+1;
        end
    end
    Lhatir(r)=Lhatir(r)/norm(subsyst.Fs(r,:)'*subsyst.Xishat,Inf);
    diff_L(r)=Lhatir(r)-(norm(subsyst.Fs(r,:)',Inf)*delta/norm(subsyst.Fs(r,:)'*subsyst.Xishat,Inf));
end
Lhat=min(diff_L);
%%%%%%%%%






    c = [-Lhat+eps;beta-1+eps;alfa-1+eps];
    ceq = [];
end
