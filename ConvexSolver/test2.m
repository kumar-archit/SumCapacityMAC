pyx = [0.4,0.1,0.5;0.3,0.2,0.5;0.5,0.4,0.1;0.2,0.799,0.001];
x0 = [0.1,0.1,0,0,0,0];
A1=eye(6);
A=-1*A1;
b = [0,0,0,0,0,0];
Aeq = [0,0,1,1,1,1;0,0,1,1,0,0;0,0,0,0,1,1;0,0,1,0,1,0;0,0,0,1,0,1;0,0,0,0,0,0];
lb = [0,0,0,0,0,0];
ub = [1,1,1,1,1,1];
 i=0.3;
 j=0.6;
for a=0.01:0.1:0.99
     for i=0.01:0.1:0.99
          for j=0.01:0.1:0.99
            p1 = [i,1-i];
            p2 = [j,1-j];
            beq = [1,p1(1),p1(2),p2(1),p2(2),0];
            fun = @(x)prog(x,a);
            nonlcon = @(x)mycon(x,pyx,p1,p2);
            options = optimoptions('fmincon','Display','off');
            % options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
             x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
             fprintf('%f,%f\n', x(1),x(2))
          end 
      end    
end    
