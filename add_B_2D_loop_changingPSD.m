function [ff_length,ini_fil,a_points,line,aux_fil] = add_B_2D_loop_changingPSD(t,ff_length,S,P)
%     create random points in a box containing the spine and select the
%     points that are inside the spine
    aux = abs(P.h_PSD-P.h_neck)/2-P.delta_S;
    aux = -aux + aux*2*rand(1000,1);
    aux2 = max(max(S(:,1)),abs(min(S(:,1))))-P.delta_S;
    aux2 = -aux2 + aux2*2*rand(1000,1);
    aux_in_all =  [aux2 aux];
    
    in = inpolygon(aux_in_all(:,1),aux_in_all(:,2),S(:,1),S(:,2));
    aux_in_all = aux_in_all(in>0,:);
%  calculate the distance to the PSD and asign the corresponding 
% probability to the points according to this distance
    d = abs((S(P.index2(2),2)-S(P.index2(1),2)).*aux_in_all(:,1)...
        -(S(P.index2(2),1)-S(P.index2(1),1)).*aux_in_all(:,2)...
        +S(P.index2(2),1)*S(P.index2(1),2)-S(P.index2(2),2)*S(P.index2(1),1))...
        ./sqrt((S(P.index2(2),2)-S(P.index2(1),2))^2+(S(P.index2(2),1)-S(P.index2(1),1))^2);  
    aux_in_all = aux_in_all(d>=P.d_arp23,:);
    d = d(d>=P.d_arp23);
    aux_d = table(d,aux_in_all);
    aux_d = unique(aux_d);
    aux_d = aux_d{:,:};
    aux_exp = exp(-aux_d(:,1)/P.tau_prob);
    aux_exp2 = ones(size(aux_exp))'*aux_exp;
    aux_d = [aux_exp/aux_exp2 aux_d];    

   aux_ind2 = [];
    while isempty(aux_ind2)
        aux_in = rand;
        aux_cumulative = cumsum(aux_d(:,1));
        aux_in = find(aux_cumulative >= aux_in,1);
        if isempty(aux_in)
            aux_in = length(aux_cumulative);
        end
        aux_ind = setdiff(P.index,P.ini_fil);
        aux_d_S = sqrt((S(aux_ind,1)-aux_d(aux_in,3)).^2 + (S(aux_ind,2)-aux_d(aux_in,4)).^2);
        aux_ind2 = find(aux_d_S <= P.d_mem);
    end
%     select the polymerization direction towards the spine membrane
    aux_ini_fil = randi(length(aux_ind2));
    aux_ini_fil = aux_ind(aux_ind2(aux_ini_fil));
    ini_fil = [P.ini_fil;aux_ini_fil];
    aux_l = aux_d(aux_in,3:4);
    a_points = [P.a_points;aux_l];
    m = (S(aux_ini_fil,2)-aux_l(2))/(S(aux_ini_fil,1)-aux_l(1));
    if S(aux_ini_fil,1)>=0
       cuad = 1;
    else
       cuad = 2;
    end
    line = [P.line;m -1 -m*aux_l(1)+aux_l(2) cuad];
    
    aux_fil = [aux_ini_fil 0 1];
    ff_length = [ff_length; aux_l t];

end
