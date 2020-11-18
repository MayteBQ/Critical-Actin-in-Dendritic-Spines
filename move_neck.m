function [S,P] = move_neck(S,P)
% check is the neck gets wider or narrower
    aux_rand = rand(2,1) < P.move_neck*P.delta_t;
    if sum(aux_rand) > 0
        if P.index3(2) < size(S,1)
            if P.index3(2)-P.index3(1) < 5
                P.index3(1) = P.index3(1) - aux_rand(2);
                P.index3(2) = P.index3(2) + aux_rand(2);
            else
                P.index3(1) = P.index3(1) + aux_rand(1) - aux_rand(2);
                P.index3(2) = P.index3(2) - aux_rand(1) + aux_rand(2);
            end   
            P.index = [1:(P.index2(1)-1) (P.index2(2)+1):(P.index3(1)-1) (P.index3(2)+1):size(S,1)]';
        end
        P.index_n = (P.index3(1)+1):(P.index3(2)-1); 
    end
    S(P.index3(1):P.index3(2),2) = P.h_neck;    
end
    