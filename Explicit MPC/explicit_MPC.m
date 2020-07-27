A = [ 0.5403, -0.8415; 0.8415,  0.5403];
B = [ -0.4597; 0.8415];
C = [1 0];
D = 0;
sys = ss(A,B,C,D,1);
model = LTISystem(sys);
model.y.min = -10;
model.y.max = 10;
model.u.min = -1;
model.u.max = 1;
model.x.penalty = QuadFunction(eye(2));
model.u.penalty = QuadFunction(1);
Tset = model.LQRSet;
PN = model.LQRPenalty;
model.x.with('terminalSet');
model.x.terminalSet = Tset;
model.x.with('terminalPenalty');
model.x.terminalPenalty = PN;
ctrl = MPCController(model,5);
expmpc = ctrl.toExplicit();
figure(1)
expmpc.partition.plot
figure(2)
expmpc.feedback.fplot
figure(3)
expmpc.cost.fplot

x0=[-2.295;-3.344];
x=x0;
for i=1:100
    tic
    uonline = ctrl.evaluate(x0);
    timeonline=toc;
    
    tic
    [log,reg]=expmpc.partition.contains(x(:,end));
    if log==1
        law=expmpc.feedback.Set(reg(1)).Functions('primal');
        useq=law.F*x(:,end)+law.g;
        u=useq(1);
    else
        break;
    end
    timeexpl=toc;
    disp(sprintf('Online time %f explicit time %f',timeonline,timeexpl));
    
%    [u, feasible, openloop] = expmpc.evaluate(x(:,end));
   x=[x A*x(:,end)+B*u];
end
figure(1)
hold on
plot(x(1,:),x(2,:),'b','LineWidth',4)