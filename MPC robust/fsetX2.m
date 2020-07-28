function [temp HH KK set_X]=fsetX2(Ad,Bd,K,Ured,Xred,xs,us)
    % -ubounds<=-Kx<=ubounds -xbounds<=x<=xbounds
    % UredH*u<=UredK
    % UredH*(-K(x-xs)+us)<=UredK
    % UredH*-K*x<=UredK-(UredH*K*xs+UredH*us)
    % XredH*x<=XredK
    
    HKtemp=Ured.H;
    UredH=HKtemp(:,1:size(Bd,2));
    UredK=HKtemp(:,end);
    HKtemp=Xred.H;
    XredH=HKtemp(:,1:size(Ad,2));
    XredK=HKtemp(:,end);
    setXU=Polyhedron('A',[UredH*-K;XredH],'b',[UredK-(UredH*K*xs+UredH*us);XredK]);
%     setXU=Polyhedron('A',[-K;K;eye(2);-eye(2)],'b',[ubounds-K*xs-us;ubounds+K*xs+us;xbounds;xbounds;xbounds;xbounds]);
    set_X{1}=setXU;

    hold on
    for i=1:100
        if i>1 && abs(volume(set_X{i-1})-volume(set_X{i}))<1e-3
            break
        end
%         volume(set_X{i})
        try
            set_X{i+1}=inv(Ad-Bd*K)*(set_X{i}-Bd*K*xs-Bd*us)&set_X{i};
    %             plot(set_X{i+1},'Color','r');
               % keyboard
        catch
            break
        end
    end
    %keyboard

    temp=set_X{end}.H;
    HH=temp(:,1:end-1);
    KK=temp(:,end);
end