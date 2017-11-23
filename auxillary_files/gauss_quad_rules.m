function [quadpts,wts] = gauss_quad_rules(n)
 
if ( n <= 20 )
    load quadtable;
%     fprintf('n = %i',n);
    quadpts = quadtable.pts{n};
    wts = quadtable.wts{n};
else
    syms x
    % calculating the quadrature points
    m = n-1;
    P0 = 1;
    P1 = x;
    for i = 1:m  % i=1 ==> P_2, so, i=n-1 ==> P_{n-1+1} = P_n
        Pn = ((2*i+1)*x*P1 - i*P0)/(i+1);  % the recurrence relation
        P0 = P1;
        P1 = Pn;
    end
    if n==1
        Pn = P1;
    end
    Pn = expand(Pn);
    quadpts=solve(vpa(Pn,32));
    quadpts=sort(quadpts);
    
    % calculating the weights
    P0 = 1;
    P1 = x;
    m = n;  %each w_i needs P_{n+1}, so we go up till P_{n+1}
    for i = 1:m  % i = 1 ==> P_2, so, i = ? ==> P_{n+1} (i = n)
        Pn = ((2*i+1)*x*P1 - i*P0)/(i+1);  % the recurrence relation
        P0 = P1;
        P1 = Pn;
    end
    wts = zeros(1,n);
    for k = 1:n
        wts(k) = vpa(2*(1-quadpts(k)^2)/(n+1)^2/subs(Pn,x,quadpts(k))^2,32);
    end
    
end

