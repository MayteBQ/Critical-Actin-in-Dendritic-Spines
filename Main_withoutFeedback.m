%%%%%%% 2D simulation with changes in PSD and neck without membrane feedback in the
%%%%%%% barbed end branching rate
rng(2012)
load('spine_steady.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% some parameters have different names in the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P.cte= phi, proportionality constant 
% P.tau_n = nucleation rate
% P.tau_prob = lambda, nucleation distance parameter
% P.tau_1 = 1/gamma_uncap, uncapping rate for - ends
% P.tau_2 = 1/gamma_sever, severing rate for - ends
% P.gamma =  gamma_cap, cappin rate for + ends

% new parameters for the branching rate 
P.FF = 7;%constant membrane force
P.BB = 4;%constant number of barbed ends
P.F_B_cte = exp(-P.FF*P.delta_y/(P.k_bT*P.BB))/P.BB;%constant branching rate

aux_t = P.t_initial:P.delta_t:P.t_end;

% for storing the data
aux_S = cell(length(aux_t),1);%mesh
aux_area = zeros(length(aux_t),1);%area
aux_foci = zeros(length(aux_t),1);%number of focus
aux_B = zeros(length(aux_t),1);%total number of barbed ends
aux_index2 = zeros(length(aux_t),2);%vertex index corresponding to the start and end of the PSD
aux_index3 = zeros(length(aux_t),2);%vertex index corresponding to the start and end of the neck

% initiate simulation, select the initial polymerization foci
[P.ini_fil,P.a_points,P.line] = initial_B_2D(S,P);
% select number of barbed ends at each focus and the state of the plus and
% minus end. 0=uncapped, 1=capped
filaments = [];
B = randi(P.B_0,length(P.ini_fil),1);
for l = 1:length(P.ini_fil)
    filaments = [filaments;P.ini_fil(l)*ones(B(l),1) zeros(B(l),1) ones(B(l),1)];
end
filaments = cap_ends(filaments,P);

ff_length = [P.a_points aux_t(1)*ones(size(P.ini_fil))];
aux_foci(1) = length(P.ini_fil);
aux_area(1) =area_S(S);
aux_S{1} = S;
aux_fil = store_data(filaments,P);%aid to trace number of barbed ends per foci
aux_B(1) = sum(aux_fil);
aux_index2(1,:) = P.index2;
aux_index3(1,:) = P.index3;
    

for j=2:length(aux_t)
               
    [ff_length,P.ini_fil,P.a_points,P.line,filaments] = nucleation_loop_2D_changingPSD(aux_t(j),ff_length,filaments,S,P);
% branch barbed ends with a constant rate
    [aux_new] = create_filaments_cte(filaments,P);
    
    [filaments] = remove_old(filaments,P); 
    [filaments] = cap_ends(filaments,P);
    filaments = [filaments; aux_new];
        
    [S,P] = move_PSD(S,P);
    [S,P] = move_neck(S,P);
        
    aux_n = unique(filaments(:,1));
    aux_dif = 1:length(P.ini_fil);
    if ~isempty(aux_n) 
        aux_nn = [];
        for ll = 1:length(aux_n)
            aux_nn = [aux_nn;find(P.ini_fil == aux_n(ll))];
        end
        C = setdiff(aux_dif,aux_nn)';
        ff_length = [ff_length;P.a_points(C,:) aux_t(j)*ones(size(C))];
        P.a_points = P.a_points(aux_nn,:);
        P.ini_fil = P.ini_fil(aux_nn);
        P.line = P.line(aux_nn,:);

    else
        ff_length = [ff_length;P.a_points aux_t(j)*ones(length(P.ini_fil),1)];
        P.a_points = [];
        P.ini_fil = [];
        P.line = [];

    end
        
    S = solve_system_threshold_movingPSD_RK(S,aux_fil,1,P);
    [S_new,P] = remesh_2D(S,P);
        
    P.K = length(S_new);
    [P.ini_fil,filaments] = move_index_2D(S_new,filaments,P);
    S = S_new;
    aux_fil = store_data(filaments,P);
        
    aux_area(j) = area_S(S);
    aux_foci(j) = length(P.ini_fil);
    aux_B(j) = sum(aux_fil);
    aux_S{j,1} = S;
    aux_index2(j,:) = P.index2;
    aux_index3(j,:) = P.index3;  
end

