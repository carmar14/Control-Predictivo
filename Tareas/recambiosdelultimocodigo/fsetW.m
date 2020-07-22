% function [temp HH KK temp2 HH2 KK2 set_w]=fsetW(K,lambda,L,ubounds,xbounds,Aw)
function [temp HH KK temp2 HH2 KK2 set_w]=fsetW(K,lambda,L,ubounds,xbounds,Aw,Mtetha)
% -ubounds<=-Kx<=ubounds -xbounds<=x<=xbounds
%     setXU=Polyhedron('A',[-K;K;eye(2);-eye(2)],'b',[ubounds-K*xs-us;ubounds+K*xs+us;xbounds;xbounds;xbounds;xbounds]);

%Calculando W_lambda --------------------------------------------
% setXt=Polyhedron('A',[eye(2) zeros(2,2);-eye(2) zeros(2,2);zeros(2,2)...
%     lambda*eye(2);zeros(2,2) -lambda*eye(2);K L;-K -L],'b',...
%     [xbounds;xbounds;xbounds;xbounds;...
%     lambda*xbounds;lambda*xbounds;lambda*xbounds;lambda*xbounds;...
%     ubounds;ubounds;ubounds;ubounds]);
 
setXt=Polyhedron('A',[eye(2) zeros(2,2);-eye(2) zeros(2,2);zeros(4,2)...
    Mtetha;zeros(4,2) -Mtetha;K L;-K -L],'b',...
    [xbounds;xbounds;xbounds;xbounds;...
    lambda*xbounds;lambda*xbounds;lambda*ubounds;lambda*ubounds;...
    
    lambda*xbounds;lambda*xbounds;lambda*ubounds;lambda*ubounds;...
    
    ubounds;ubounds;ubounds;ubounds]);


%%% Comentarios: aqui falta algo. Mirando a la definicion de Wlambda, se
%%% pide que Mtheta*Theta pertenezca a lambda Z. Mtheta*Theta incluye xs y
%%% us. Entonces no es correcto lo que escribiste. Hay que poner
%%% lambda*eye(2)*Mthetax y -lambda*eye(2)*Mthetax y tambien
%%% lambda*eye(2)*Mthetau y -lambda*eye(2)*Mthetau con Mthetax las
%%% primeras 2 lineas de Mtheta y Mthetau las siguentes

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

%     figure
%     plot(set_w{end}.projection([1 2]))
%     hold on

temp2=set_w{end}.projection([1 2]).H;
HH2=temp2(:,1:end-1);
KK2=temp2(:,end);

%%% Comentarios: quidado! Las restricciones del problema de optimizacion de
%%% la seccion 3 necesitan el intero conjunto y no solo su proyeccion en
%%% X!!

    temp=set_w{end}.H;
    HH=temp(:,1:end-1);
    KK=temp(:,end);
end