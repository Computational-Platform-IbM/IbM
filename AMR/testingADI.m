a = 1; b = 2; c = 3;

A = [b c 0 0 0 0 a;...
a b c 0 0 0 0;...
0 a b c 0 0 0;...
0 0 a b c 0 0;...
0 0 0 a b c 0;...
0 0 0 0 a b c;...
c 0 0 0 0 a b;];

d = rand(7,1);

n = 1000000;

tic
for i = 1:n
    x_SM = ShermanMorrison(A, d);
end
toc

tic
for i = 1:n
    x_OG = OG(A, d);
end
toc

tic
for i = 1:n
    x_ST = Stavanger(A, d);
end
toc

[x_OG x_SM x_ST]


function x = ShermanMorrison(A, d)
    u = [-A(1,1) 0 0 0 0 0 A(end, 1)]';
    v = [1 0 0 0 0 0 -A(1, end)/A(1,1)]';
    A_star = A - u*v';
    y = A_star \ d;
    q = A_star \ u;
    
    x = y - ((v' * y)/(1 + v' * q)) * q;
end

function x = OG(A, d)
    x = A \ d;
end

function x = Stavanger(A, d)
    rhs = d(1:end-1);
    Ac = A(1:end-1, 1:end-1);
    x1 = Ac \ rhs;
    q = zeros(length(rhs), 1);
    q(1) = A(1, end); q(end) = A(end - 1, end);
    x2 = Ac \ -q;
    
    
    x_end = (d(end) - A(end, 1)*x1(1) - A(end, end-1)*x1(end)) / (A(end, end) + A(end, 1)*x2(1) + A(end, end-1)*x2(end));
    x = [x1 + x2*x_end; x_end];
end