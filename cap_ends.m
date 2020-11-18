function [filaments] = cap_ends(filaments,P)
%     capping of the plus end
    aux = find(filaments(:,2)==0);
    aux_rand = rand(size(aux)) < P.gamma*P.delta_t;
    %     if the plus end in capped then the filament is neglected 
    filaments(aux(aux_rand),:) = [];  

    %     uncapping the minus end 
    aux= find(filaments(:,3)==1);
    aux_rand = rand(size(aux)) < P.delta_t/P.tau_1;
    filaments(aux(aux_rand),3) = 0;
end