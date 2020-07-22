function [temp HH KK set_w]=fsetW(K,lambda,L,ubounds,xbounds,Aw)
% -ubounds<=-Kx<=ubounds -xbounds<=x<=xbounds
%     setXU=Polyhedron('A',[-K;K;eye(2);-eye(2)],'b',[ubounds-K*xs-us;ubounds+K*xs+us;xbounds;xbounds;xbounds;xbounds]);

%Calculando W_lambda --------------------------------------------
setXt=Polyhedron('A',[eye(2) zeros(2,2);-eye(2) zeros(2,2);zeros(2,2)...
    lambda*eye(2);zeros(2,2) -lambda*eye(2);K L;-K -L],'b',...
    [xbounds;xbounds;xbounds;xbounds;...
    lambda*xbounds;lambda*xbounds;lambda*xbounds;lambda*xbounds;...
    ubounds;ubounds;ubounds;ubounds]);

%------------------------------------------------

set_w{1}=setXt;
volume(set_w{1})

for i=1:100
    if i>1 && abs(volume(set_w{i-1}.projection([1 2]))-volume(set_w{i}.projection([1 2])))<1e-3
        
        break
    end
    volume(set_w{i}.projection([1 2]))
    try
        set_w{i+1}=inv(Aw)*(set_w{i})&set_w{i};%-Bd*K*xs-Bd*us)&set_X{i};
        
        %             plot(set_X{i+1},'Color','r');
        % keyboard
    catch
        
        disp("try")
        break
        
    end
end
%keyboard
i
%     figure
%     plot(set_w{end}.projection([1 2]))
%     hold on
temp=set_w{end}.projection([1 2]).H;
HH=temp(:,1:end-1);
KK=temp(:,end);
%     temp=set_w{end}.H;
%     HH=temp(:,1:end-1);
%     KK=temp(:,end);
end