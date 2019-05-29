%mycon.m

function [c,ceq]=mycon(x,pyx,p1,p2)
func1 = @(x)prog1(x,p2,pyx);
func2 = @(x)prog2(x,p1,pyx);
func12 = @(x)prog12(x,pyx);
c=[x(1) - func1(x); x(2) - func2(x); x(1)+x(2) - func12(x)];
ceq=[];