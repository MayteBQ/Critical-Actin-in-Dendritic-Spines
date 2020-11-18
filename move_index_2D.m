function [aux_ini_fil,filaments] = move_index_2D(S_new,filaments,P)
% update the vertex corresponding to the polymerization direction 
    aux_ini_fil = zeros(size(P.ini_fil));    
    for l = 1:length(P.ini_fil)
        if P.line(l,4) == 1
            ind = find(S_new(P.index,1) >= 0);
            aux_d = abs(P.line(l,1)*S_new(P.index(ind),1) + P.line(l,2)*S_new(P.index(ind),2) + P.line(l,3))./...
                sqrt(P.line(l,1)^2 + P.line(l,2)^2);
        else
            ind = find(S_new(P.index,1) < 0);
            aux_d = abs(P.line(l,1)*S_new(P.index(ind),1) + P.line(l,2)*S_new(P.index(ind),2) + P.line(l,3))./...
                sqrt(P.line(l,1)^2 + P.line(l,2)^2);
        end
        aux_ini_fil(l) = P.index(ind(aux_d == min(aux_d)));
        filaments(filaments(:,1)==P.ini_fil(l),1) = aux_ini_fil(l);
    end
end