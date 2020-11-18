function [aux_new] = create_filaments(filaments,a,B,P)
%     creation of a new filament at the membrane if + end uncapped
    aux = find(filaments(:,2)==0);
    aux_rand = rand(size(aux)) < P.delta_t*P.cte*P.delta_y*P.k_on*(a.*exp(-P.w(filaments(aux,1))./(B(filaments(aux,1)))))./B(filaments(aux,1));
    aux_new= [filaments(aux(aux_rand),1), zeros(size(aux(aux_rand))), ....
        ones(size(aux(aux_rand)))];
    if ~isempty(aux_new)
        [aux_new] = cap_ends(aux_new,P);
    end
end