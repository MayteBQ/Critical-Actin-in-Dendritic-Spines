function S_int = solve_system_threshold_movingPSD_RK(S,aux_fil,actin_dyn,P)
    delta_t = P.delta_t;
    count = 0;
% Solve the system, and fix the vertices corresponding to the PSD and neck
    S_int = solve_system_RK(S,aux_fil,actin_dyn,delta_t,P);
    S_int(P.index3,:) = S(P.index3,:);
    S_int(P.index_n,:) = S(P.index_n,:);
    S_int(P.index2,:) = S(P.index2,:);
    S_int(P.index_psd,:) = S(P.index_psd,:);
%     Check how much the vertices move, and if it is more than a threshold
%     then slide the time integration windos
    d_max = max(sqrt((S_int(:,1)-S(:,1)).^2 + (S_int(:,2)-S(:,2)).^2));
    while d_max >  P.d_max_int
        delta_t = delta_t/2;
        S_int = solve_system_RK(S,aux_fil,actin_dyn,delta_t,P);
        if (P.index3(2)+2) <= P.index(end) 
            S_int((P.index3(1)-2):(P.index3(2)+2),:) = S((P.index3(1)-2):(P.index3(2)+2),:);
        elseif (P.index3(2)+1) == P.index(end) 
            S_int((P.index3(1)-2):(P.index3(2)+1),:) = S((P.index3(1)-2):(P.index3(2)+1),:);
            S_int(1,:) = S(1,:);
        elseif P.index3(2) == P.index(end) 
            S_int((P.index3(1)-2):(P.index3(2)),:) = S((P.index3(1)-2):(P.index3(2)),:);
            S_int(1:2,:) = S(1:2,:);
        end
        S_int(P.index2,:) = S(P.index2,:);
        S_int(P.index_psd,:) =  S(P.index_psd,:);
        d_max = max(sqrt((S_int(:,1)-S(:,1)).^2 + (S_int(:,2)-S(:,2)).^2));
        count = count+1;

    end

    if count > 0
        for j=1:(2^count-1)
            S_int_new = solve_system_RK(S_int,aux_fil,actin_dyn,delta_t,P);
            if (P.index3(2)+2) <= P.index(end)
                S_int_new((P.index3(1)-2):(P.index3(2)+2),:) = S((P.index3(1)-2):(P.index3(2)+2),:);
            elseif (P.index3(2)+1) == P.index(end) 
                S_int_new((P.index3(1)-2):(P.index3(2)+1),:) = S((P.index3(1)-2):(P.index3(2)+1),:);
                S_int_new(1,:) = S(1,:);
            elseif P.index3(2) == P.index(end) 
                S_int_new((P.index3(1)-2):(P.index3(2)),:) = S((P.index3(1)-2):(P.index3(2)),:);
                S_int_new(1:2,:) = S(1:2,:);
            end
            S_int_new(P.index2,:) = S_int(P.index2,:);
            S_int_new(P.index_psd,:) =S_int(P.index_psd,:);
            S_int = S_int_new;
        end
    end
   
end


    
    
        
