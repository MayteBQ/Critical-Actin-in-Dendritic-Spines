function [S,P] = move_PSD(S,P)
    aux_rand = rand(4,1) < P.move_PSD*P.delta_t;
%     check if the PSD moves to the left or to the right
    if sum(aux_rand) > 0
        P.index2(1) = P.index2(1) - aux_rand(1) + aux_rand(2);
        P.index2(2) = P.index2(2) + aux_rand(3) - aux_rand(4);
        if P.index3(2) <= size(S,1)
            P.index = [1:(P.index2(1)-1) (P.index2(2)+1):(P.index3(1)-1) (P.index3(2)+1):size(S,1)]';
        else
            P.index = [1:(P.index2(1)-1) (P.index2(2)+1):(P.index3(1)-1)]';
        end
        P.index_psd = (P.index2(1)+1):(P.index2(2)-1);
    end
%     check if the PSD moves up or down
    aux_rand = rand(2,1) < P.move_PSD_h*P.delta_t;
    if sum(aux_rand) > 0
        P.h_PSD = P.h_PSD + aux_rand(1)*P.d_move_PSD - aux_rand(2)*P.d_move_PSD;
    end
    S(P.index2(1):P.index2(2),2) = P.h_PSD;    
end
    
