clear;
clc;
f = @(x,y) x-(1/y);

y=1;
yn=1;
h=0.1;
for x = [0 0.1]
    j=0;
    yi = y + h*f(x,y);
    while( abs(yi-y)/abs(yi) > 0.0001)
        y=yi;
        yi = yn + h*(f(x,yn)+ f(x+h,y))/2;
        j=j+1;
    end
    fprintf("Inner Function converge after %d steps\n",j);
    yn=yi;
end