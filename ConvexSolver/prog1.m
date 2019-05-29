%prog1.m
function f1=prog1(x,p2,pyx)
f1=0.0;
for y = 1:3
    for i = 1:2
        for j = 1:2
            sum1=0.0;
            for k=1:2
                sum1 = sum1 + x(2*k+j)*pyx(2*(k-1)+j,y);
            end
          f1 = f1 + x(2*i+j)*pyx(2*(i-1)+j,y)*log(pyx(2*(i-1)+j,y)*p2(j)/sum1);  
        end
    end
end
