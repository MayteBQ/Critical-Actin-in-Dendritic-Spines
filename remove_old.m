function [filaments] = remove_old(filaments,P)
      %     remove the old uncapped minus end 
    aux = find(filaments(:,3)==0); 
    aux_rand = rand(size(aux)) < P.delta_t/P.tau_2;
    filaments(aux(aux_rand),:) = [];
end