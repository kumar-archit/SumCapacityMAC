function f2=prog2(x,p1,pyx)
f2=0.0;
for y = 1:3
    for i = 1:2
        for j = 1:2
            sum2=0.0;
            for k=1:2
                sum2 = sum2 + x(2*i+k)*pyx(2*(i-1)+k,y);
            end
          f2 = f2 + x(2*i+j)*pyx(2*(i-1)+j,y)*log(pyx(2*(i-1)+j,y)*p1(i)/sum2);  
        end
    end
end
