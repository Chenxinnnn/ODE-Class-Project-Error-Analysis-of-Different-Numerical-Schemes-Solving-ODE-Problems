sg = @(y) 1/(1+exp(-y));
rk = rk4(0.5,0.1,100,0);
t = rk(2);
e = sg(t) - rk(1);

x=linspace(0,t,100);
plot(x,sg(x))
plot(x,rk(3))

function [y_list,t] = rk4(y0,h,n,t0)
    y = y0;
    y_list = zeros(n,1);
    for i = 1:n
        L_1 = f(y);
        L_2 = f(y + 0.5*h*L_1);
        L_3 = f(y + 0.5*h*L_2);
        L_4 = f(y + h*L_3);
        y = y + h/6*(L_1 + 2*L_2 + 2*L_3 + L_4);
        y_list(i) = y;
        t = t0 + n*h;
    end
end
    
function f = f(x)
    f = x-x^2;
end

