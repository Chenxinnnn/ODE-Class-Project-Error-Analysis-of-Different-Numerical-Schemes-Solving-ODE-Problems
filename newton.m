clear all
close all
f=@(y) y-y^2;
df=@(y) 1-2*y; %derivative of f (in a systems of ODE, it will be 
%the jacobian matrix of f)


ts=100; %time steps


k=20/ts;
y0=0.5; %initial condition

%function to set to zero:
%y will be the solution at time n+1   
%yn is the current y
F=@(y,yn) y-yn-k*f(y);
J=@(y) 1 - k*df(y);%jacobian matrix for F. In this case is just a number

y=NaN(1,ts+1);%initialize the solution to NaN
y(1,1)=y0;
tol=k^2; %using backward Euler

for n=1:ts
    y(1,n+1)=y(1,n); %initial guess for Newton's method

    res=-J(y(1,n+1))\F(y(1,n+1),y(1,n)); %Deltax
    while (norm(res,inf)>tol)
        y(1,n+1)=y(1,n+1)+res;
        res=-J(y(1,n+1))\F(y(1,n+1),y(1,n));
    end
    y(1,n+1)=y(1,n+1)+res;
 end

t=linspace(0,20,ts+1);
sgg = zeros(1,101);
e = zeros(1,101);
for i=1:100
    sgg(i) = sg(t(i));
    e(i) = abs(y(i)-sgg(i));
end
plot(t,y,t,sgg,t,e)
xlabel('t')
ylabel('y')
legend({'Newton Method','sigmoid function','error'})
axis([0,10.2,-0.1,1.2])

function f = sg(x)
    f = 1 / (1 + exp(-x));
end