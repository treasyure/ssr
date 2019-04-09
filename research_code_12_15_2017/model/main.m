load('matlab.mat')
[t,y]=ode23s(@(t,y)odeEquation(t,y,a),[0 6*30+40],[FSH LH FSHp LHp phi omega lamda S Ty T E2 P4]);
%%plot(t,y(:,3))
axis([0 100 0 300])
init=y(end,:);
y2= [155.1824    4.2012   84.6065  207.0007    0.3776    9.8888    5.6286    0.0432    0.0031  259.5425   41.2947    1.1785];