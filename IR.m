function [f,S,D,O] = IR(flag,A,B,P)
% Calculate the shape descriptors
    aux =1:P.n-1;
    f = zeros(P.n*2,1);
    S = [];
    D = [];
    O = [];
    if flag == 1
        for j = 1:P.n*2
            f(j) = A(1) + sum(A(2:P.n).*cos(P.theta(j).*aux) + ...
                B(2:P.n).*sin(P.theta(j).*aux)) + A(P.n+1)*cos(P.n*P.theta(j));
        end
    else
        S = A(1);
        D = zeros(P.n*2,1);
        O = zeros(P.n*2,1);
        for j = 1:P.n*2
            f(j) = A(1) + sum(A(2:P.n).*cos(P.theta(j).*aux) + ...
                B(2:P.n).*sin(P.theta(j).*aux)) + A(P.n+1)*cos(P.n*P.theta(j));
            D(j) = sum(A(2:2:P.n).*cos(P.theta(j).*aux(1:2:end)) + ...
                B(2:2:P.n).*sin(P.theta(j).*aux(1:2:end))) + A(P.n+1)*cos(P.n*P.theta(j));
            O(j) = sum(A(3:2:P.n).*cos(P.theta(j).*aux(2:2:end)) + ...
                B(3:2:P.n).*sin(P.theta(j).*aux(2:2:end)));
        end
    end
    
end