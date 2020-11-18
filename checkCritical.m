% cell array with different simulations
% each data set contains P aux_B aux_S aux_area aux_foci aux_index2
% aux_index3 aux_ff_length

Data_spines = {'Sim1.mat','Sim2.mat','Sim3.mat','Sim4.mat','Sim5.mat',...
    'Sim6.mat','Sim7.mat','Sim8.mat','Sim9.mat','Sim10.mat'};


avalanche_all = [];
slope_s = zeros(length(Data_spines),1);
int_T2 = [1 5 10 20 35];
avalanche_T2_all = cell(size(int_T2));
avalanche_rest_all = [];
avalanche_rest2_all = [];
data_evo.P.delta_t = 1/8;
T = data_evo.P.delta_t:data_evo.P.delta_t:2000*data_evo.P.delta_t;
aux_F_all = cell(size(T,1),length(Data_spines));
sigma_all = cell(length(Data_spines),1);
t_sigma_all = cell(length(Data_spines),1);


for ss = 1:length(Data_spines)
    data_evo = load(Data_spines{1,ss});
    data_evo.P.t_initial = 10*60;
    aux_t = data_evo.P.t_initial:data_evo.P.delta_t:data_evo.P.t_end;
    ini_t = data_evo.P.t_initial/data_evo.P.delta_t + 1;
    %%%%%%%%%%%%%%
    %size power-law distribution 
    %%%%%%%%%%%%%
    avalanche = [];
    sigma = [];
    t_sigma = [];
    B_ends = data_evo.aux_B(ini_t:end);
    j = 1;
    while j <= length(aux_t)
        if B_ends(j) > 0
            aux_a = B_ends(j);
            j_0 = j;
            j = j + 1;
            while j < length(aux_t) && B_ends(j) > 0 
                aux_a = aux_a + B_ends(j);
                sigma = [sigma;B_ends(j)/B_ends(j-1)];
                t_sigma = [t_sigma;aux_t(j)];
                j = j + 1;
            end
%             save avalanche data: number of barbed end, furation in time
%             bins, inital and final time
            avalanche = [avalanche; aux_a (j-1-j_0) aux_t(j_0) aux_t(j-1)];
        else 
            j = j + 1;
        end
    end
    avalanche = avalanche(2:end-1,:);
    avalanche_all = [avalanche_all;avalanche];
    sigma_all{ss,1} = sigma;
    t_sigma_all{ss,1} = t_sigma;


    %%%%%%%%%%%%%%
    %temporally scale free
    %%%%%%%%%%%%%

    for kk = 1:length(int_T2)
        avalanche_T2 = [];
        T2 = data_evo.P.t_initial:data_evo.P.delta_t*int_T2(kk):data_evo.P.t_end;
        B_ends_T2 = zeros(size(T2));
        B_ends_T2(1) = B_ends(1);
        for j = 2:length(T2)
            ini = find(aux_t > T2(j-1),1);
            fin = find(aux_t == T2(j));
            B_ends_T2(j) = sum(B_ends(ini:fin));
        end

        j = 1;
        while j <= length(T2)
            if B_ends_T2(j) > 0
                aux_a = B_ends_T2(j);
                j_0 = j;
                j = j + 1;
                while j <= length(T2) &&  B_ends_T2(j) > 0
                    aux_a = aux_a + B_ends_T2(j);
                    j = j + 1;
                end
                avalanche_T2 = [avalanche_T2; aux_a (j-1-j_0) T2(j_0) T2(j-1)];
            else 
                j = j + 1;
            end
        end
        avalanche_T2 = avalanche_T2(2:end-1,:);
        avalanche_T2_all{1,kk}=[avalanche_T2_all{1,kk};avalanche_T2];

    end


    %%%%%%%%%%%%%%
    %spatially scale-free 
    %%%%%%%%%%%%%
   
    B_mean = B_ends./ data_evo.aux_foci(ini_t:end);
    B_mean( data_evo.aux_foci(ini_t:end) == 0) = 0;
    aux = data_evo.ff_length;
    time = [];
    while ~isempty(aux)
        A = (aux(:,1) == aux(1,1));
        B = (aux(:,2) == aux(1,2));
        C = A.*B;
        D = aux(C==1,:);
        if sum(C) > 1
            time = [time; D(1,:) D(2,3)];
            aux = aux(C==0,:);
        elseif sum(C) == 1
            time = [time; D(1,:) aux_t(end)];
            aux = aux(C==0,:);
        end
    end
    time2 = time;
    rest_foci = zeros(length(aux_t),1);
    aux = find(time2(:,3) == aux_t(1));
    if length(aux) > 2
        aux_ind = setdiff(aux,aux(1:2));
        aux = aux(1:2);
        time2(aux_ind,:) = [];
    end
    rest_foci(1) = length(aux);
    for j = 2:length(aux_t)
        aux = find(time2(:,3) == aux_t(j));
        aux2 = find(time2(:,3) < aux_t(j) & time2(:,4)>=aux_t(j));
        if ~isempty(aux)
            if isempty(aux2)
                if length(aux) > 2
                    aux_ind = setdiff(aux,aux(1:2));
                    aux = aux(1:2);
                    time2(aux_ind,:) = [];
                end
            elseif length(aux2) == 2
                time2(aux,:) = [];
                aux = [];
            elseif length(aux2) == 1
                aux_ind = setdiff(aux,aux(1));
                aux = aux(1);
                time2(aux_ind,:) = [];
            end
        end
        rest_foci(j) = length(aux2) + length(aux);
    end
    B_rest = B_mean.*rest_foci;
        
    avalanche_rest = [];
    j = 1;
    while j <= length(aux_t)
        if B_rest(j) > 0
            aux_a = B_rest(j);
            j_0 = j;
            j = j + 1;
            while j <= length(aux_t) && B_rest(j) > 0 
                aux_a = aux_a + B_rest(j);
                j = j + 1;
            end
            avalanche_rest = [avalanche_rest; aux_a (j-1-j_0) aux_t(j_0) aux_t(j-1)];
        else
            j = j + 1;
        end
    end
    avalanche_rest = avalanche_rest(2:end-1,:);
    avalanche_rest_all = [avalanche_rest_all;avalanche_rest];

%     for different number of Foci
    time3 = time;
    rest_foci2 = zeros(length(aux_t),1);
    aux = find(time3(:,3) == aux_t(1));
    if length(aux) > 1
        aux_ind = setdiff(aux,aux(1));
        aux = aux(1);
        time3(aux_ind,:) = [];
    end
    rest_foci2(1) = length(aux);

    for j = 2:length(aux_t)
        aux = find(time3(:,3) == aux_t(j));
        aux2 = find(time3(:,3) < aux_t(j) & time3(:,4)>=aux_t(j));
        if ~isempty(aux)
            if isempty(aux2)
                if length(aux) > 1
                    aux_ind = setdiff(aux,aux(1));
                    aux = aux(1);
                    time3(aux_ind,:) = [];
                end
            elseif length(aux2) == 1
                time3(aux,:) = [];
                aux = [];
            end
        end
        rest_foci2(j) = length(aux2) + length(aux);
    end
    B_rest2 = B_mean.*rest_foci2;
    avalanche_rest2 = [];
    j = 1;
    while j <= length(aux_t)
        if B_rest2(j) > 0
            aux_a = B_rest2(j);
            j_0 = j;
            j = j + 1;
            while j <= length(aux_t) && B_rest2(j) > 0
                aux_a = aux_a + B_rest2(j);
                j = j + 1;
            end
            avalanche_rest2 = [avalanche_rest2; aux_a (j-1-j_0) aux_t(j_0) aux_t(j-1)];
        else
            j = j + 1;
        end
    end
    avalanche_rest2 = avalanche_rest2(2:end-1,:);
    avalanche_rest2_all = [avalanche_rest2_all;avalanche_rest2];

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% For Allen (Fano) Factor
    %%%%%%%%%%%%%%%%%%
    
    
    for j = 1:length(T)
        aux_T = data_evo.P.t_initial:T(j):data_evo.P.t_end;
        aux_F = zeros(length(aux_T)-1,1);
        for k = 1:length(aux_F)
            ini = find(aux_t == aux_T(k));
            fin = find(aux_t == aux_T(k+1));
            aux_F(k) = sum(B_ends(ini:fin-1));
        end
        aux_F_all{j,ss} = aux_F;
    end

end

avalanche_all_aux = avalanche_all;
avalanche_all_aux(avalanche_all_aux(:,1)==1,:) = [];
figure
[N,edges] = histcounts(avalanche_all_aux(:,1),-1+min(avalanche_all_aux(:,1)):max(avalanche_all_aux(:,1)));
prob = N./sum(N);
loglog(edges(2:end),prob,'k')
hold on 
loglog(edges(2:end),.20*edges(2:end).^(-1),'r--')
set(gca,'fontsize',14)
xlabel('size')
ylabel('P(size)')


duration = avalanche_all_aux(:,4)-avalanche_all_aux(:,3);
[N_t,edges_t] = histcounts(duration,(min(duration)-data_evo.P.delta_t):data_evo.P.delta_t:max(duration));
prob_t = N_t./sum(N_t);
figure
loglog(edges_t(2:end),prob_t,'k')
hold on 
loglog(edges_t(3:end),.025*edges_t(3:end).^(-1.36/2),'r--')
xlabel('duration')
ylabel('P(duration)')
set(gca,'fontsize',14)
xlabel('lifetime')
ylabel('P(lifetime)')

figure
s_c = 1;
aux = avalanche_all(avalanche_all(:,1)>=s_c,:);
tau = aux(2:end,3)-aux(1:end-1,4);
tau(tau < 0) = [];
[N,edges_II] = histcounts(tau,(min(tau)-data_evo.P.delta_t):data_evo.P.delta_t:max(tau));
prob_II = N./sum(N);
loglog(edges_II(2:end),prob_II,'k')
hold on 
loglog(edges_II(2:end),.04*edges_II(2:end).^(-1),'r--')
set(gca,'fontsize',14)
xlabel('\eta')
ylabel('P(\eta)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% for supplementary material
% temporal scale free
figure
for jj = 1:length(int_T2)
    aux_avalanche_T2 = avalanche_T2_all{1,jj};
    [N,edges] = histcounts(aux_avalanche_T2(:,1),0:max(aux_avalanche_T2(:,1)));
    prob = N./sum(N);
    loglog(edges(2:end),prob)
    hold on
end
loglog(edges(2:end),.1./edges(2:end),'k:','linewidth',2)
xlabel('size')
ylabel('P(size)')
set(gca,'fontsize',10)
axis('tight')
legend(['{\Delta}t = ' num2str(int_T2(1),'%2d') '{\Delta}t'],...
    ['{\Delta}t = ' num2str(int_T2(2),'%2d') '{\Delta}t'],...
    ['{\Delta}t = ' num2str(int_T2(3),'%2d') '{\Delta}t'],...
    ['{\Delta}t = ' num2str(int_T2(4),'%2d') '{\Delta}t'])

% spatially scale-free
figure;
[N,edges] = histcounts(avalanche_all_aux(:,1),0:max(avalanche_all_aux(:,1)));
prob = N./sum(N);
loglog(edges(2:end),prob)
hold on
avalanche_rest_all_aux = avalanche_rest_all;
avalanche_rest_all_aux(avalanche_rest_all_aux(:,1)==1,:) = [];
[N,edges] = histcounts(avalanche_rest_all_aux(:,1),0:max(avalanche_rest_all_aux));
prob_rest = N./sum(N);
loglog(edges(2:end),prob_rest)
[p_rest,e_rest] = polyfit(log(edges(2:end)),log(prob_rest),1);
[f_rest,err_rest] = polyval(p_rest,log(edges(2:end)),e_rest);
avalanche_rest2_all_aux = avalanche_rest2_all;
avalanche_rest2_all_aux(avalanche_rest2_all_aux(:,1)==1,:) = [];
[N,edges] = histcounts(avalanche_rest2_all(:,1),0:max(avalanche_rest2_all_aux(:,1)));
prob_rest2 = N./sum(N);
loglog(edges(2:end),prob_rest2)
[p_rest2,e_rest2] = polyfit(log(edges(2:end)),log(prob_rest2),1);
[f_rest2,err_rest2] = polyval(p_rest2,log(edges(2:end)),e_rest2);
loglog(edges(2:end),.2./edges(2:end),'k:','linewidth',2)
xlabel('size')
ylabel('P(size)')
legend('All foci','2 foci','1 focus')
set(gca,'fontsize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% inter-avalanche intervals
%%%%%%%%%%%%%%%%%%
figure;
s_c = [1 5 10 20 50];
for kk=1:length(s_c)
    aux = avalanche_all(avalanche_all(:,1)>=s_c(kk),:);
    tau = aux(2:end,3)-aux(1:end-1,4);
    tau(tau < 0) = [];
    aux_R = mean(1./tau);
    [N,edges] = histcounts(tau.*aux_R);
    prob_II = N./sum(N);
    loglog(edges(2:end),prob_II)
    [p_restII,e_restII] = polyfit(log(edges(2:end)),log(prob_II),2);
    [f_restII,err_restII] = polyval(p_restII,log(edges(2:end)),e_restII);
    hold on
end
loglog(edges(7:end),15*edges(7:end).^(-3),'k:','linewidth',2)
xlabel('s_c R(s_c)')
ylabel('P(\eta_c,s_c)/R(s_c)')
legend(['s_c = ' num2str(s_c(1),'%2d')],...
    ['s_c = ' num2str(s_c(2),'%2d')],...
    ['s_c = ' num2str(s_c(3),'%2d')],...
    ['s_c = ' num2str(s_c(4),'%2d')],...
    ['s_c = ' num2str(s_c(5),'%2d')])
set(gca,'fontsize',10)

% Allen Factor

F = zeros(size(T));
A = zeros(size(T));

for j = 1:length(T)
    FF = [];
    AA = [];
    for ss = 1:length(Data_spines)
        FF = [FF;aux_F_all{j,ss}];
        AA = [AA;(aux_F_all{j,ss}(2:end)-aux_F_all{j,ss}(1:end-1)).^2];
    end
    F(j) = var(FF)/mean(FF);
    A(j) = mean(AA)/(2*mean(FF));
end

figure;
loglog(T,A)
xlabel('T')
ylabel('Allen Factor')
hold on 
loglog(T(1:150),T(1:150).^(2),'k:','linewidth',2)
set(gca,'fontsize',10)

