function f12=prog12(x,pyx)
f12=0.0;
for y = 1:3
    for i = 1:2
        for j = 1:2
            sum12=0.0;
            for k=1:2
                for l=1:2
                    sum12 = sum12 + x(2*k+l)*pyx(2*(k-1)+l,y);
                end
            end
          f12 = f12 + x(2*i+j)*pyx(2*(i-1)+j,y)*log(pyx(2*(i-1)+j,y)/sum12);  
        end
    end
end
