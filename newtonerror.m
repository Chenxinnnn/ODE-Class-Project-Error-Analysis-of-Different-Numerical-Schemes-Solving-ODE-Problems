clear all
close all
f=@(y) y-y^2;
df=@(y) 1-2*y; %derivative of f (in a systems of ODE, it will be 
%the jacobian matrix of f)

E_max=zeros(1,10);
E_mean = zeros(1,10);
ts=100; %time steps

for r = 5:5:50
k=r/ts;
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

E_max(r/5) = max(e);
E_mean(r/5) = mean(e);
end

h = 5:5:50;
plot(h,E_mean,h,E_max)
legend({'Mean Error','Max Error'})
disp(E_mean)

function f = sg(x)
    f = 1 / (1 + exp(-x));
end