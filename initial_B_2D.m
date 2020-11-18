function [ini_fil,a_points,line] = initial_B_2D(S_0,P)

    ini_fil = [];
    a_points = [];
    line = [];
%     create random points in a box containing the inital shape and select
%     those inside the inital shape
    aux_in_all =  -P.R + P.L*rand(50000,2);
    in = inpolygon(aux_in_all(:,1),aux_in_all(:,2),S_0(:,1),S_0(:,2));
    aux_in_all = aux_in_all(in>0,:);
    
% calculate the distance to the PSD    
    d = zeros(size(aux_in_all,1),1);
    d(aux_in_all(:,1)<-P.r_PSD) = sqrt((aux_in_all(aux_in_all(:,1)<-P.r_PSD,1)+P.r_PSD).^2 +...
        (aux_in_all(aux_in_all(:,1)<-P.r_PSD,2)-P.h_PSD).^2);
    d(aux_in_all(:,1)>P.r_PSD) = sqrt((aux_in_all(aux_in_all(:,1)>P.r_PSD,1)-P.r_PSD).^2 +...
        (aux_in_all(aux_in_all(:,1)>P.r_PSD,2)-P.h_PSD).^2);
    d(aux_in_all(:,1)<=P.r_PSD & aux_in_all(:,1)>=-P.r_PSD) = ...
        abs(aux_in_all(aux_in_all(:,1)<=P.r_PSD & aux_in_all(:,1)>=-P.r_PSD,2)-P.h_PSD);

%     discard the points that are closer than d_arp23 to the PSD and
%     calculate exp(-dist/lambda)
    aux_in_all = aux_in_all(d>=P.d_arp23,:);
    d = d(d>=P.d_arp23);
    aux_d = table(d,aux_in_all);
    aux_d = unique(aux_d);
    aux_d = aux_d{:,:};
    aux_exp = exp(-aux_d(:,1)/P.tau_prob);
    aux_exp2 = ones(size(aux_exp))'*aux_exp;
    aux_d = [aux_exp/aux_exp2 aux_d];
    
%     for every focus
    while length(ini_fil) < P.ini_foci
        aux_ind2 = [];
%         select the location of the focus, here focus near the PSD are
%         more likely.
        while isempty(aux_ind2)
            aux_in = rand;
            aux_cumulative = cumsum(aux_d(:,1));
            aux_in = find(aux_cumulative >= aux_in,1);
            if isempty(aux_in)
                aux_in = length(aux_cumulative);
            end
            aux_ind = setdiff(P.index,ini_fil);
            aux_d_S = sqrt((S_0(aux_ind,1)-aux_d(aux_in,3)).^2 + (S_0(aux_ind,2)-aux_d(aux_in,4)).^2);
            aux_ind2 = find(aux_d_S <= P.d_mem); %check mesh points that are within a distance of d_mem
        end
%         select the polymerization direction and the corresponding vertex
%         index from the mesh
        aux_ini_fil = randi(length(aux_ind2));
        aux_ini_fil = aux_ind(aux_ind2(aux_ini_fil));

        if ~ismember(aux_ini_fil,ini_fil)
            ini_fil = [ini_fil;aux_ini_fil];
            aux_l = aux_d(aux_in,3:4);
            a_points = [a_points;aux_l];%location of the polymerization focus
%             slope of the line between from the polymerization focus to
%             the mesh vertex from the polymerization direction
            m = (S_0(aux_ini_fil,2)-aux_l(2))/(S_0(aux_ini_fil,1)-aux_l(1)); 
            if S_0(aux_ini_fil,1)>=0
               cuad = 1;
            else
               cuad = 2;
            end
            line = [line;m -1 -m*aux_l(1)+aux_l(2) cuad];
        end
        aux_d(aux_in,:) = []; 
    end
end