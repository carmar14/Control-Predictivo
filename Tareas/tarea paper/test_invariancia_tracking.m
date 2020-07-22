% close all
x00=[0;-2]; %[0.3; -0.1] [-1 ;-1]
x=x00;
for ii=1:100
    temp=A*x(:,ii)+B*(K*(x(:,ii)-xs(:,1))+us);
    x=[x temp];
end
figure
plot(set_w{end}.projection([1 2]),'Color','y','FaceColor','none');
hold on
plot(x(1,:),x(2,:),'*')

