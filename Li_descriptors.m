%%%%%%%%%%%%%% calculate Li descriptors
% number of sampling points
P.n = 12;
P.theta = ((1:2*P.n)-1)*pi./P.n;
P.delta_theta = P.theta(2)-P.theta(1);
P.pi_theta = find(P.theta>=pi,1);
% ray that instersects the ROI
P.r = (0.001:0.001:10)';
% for significance test
P.b_min = 0.1;
P.alpha_test = 0.1;

aux_t = P.t_initial:P.delta_t:P.t_end;

% time sampling vector
P.n_frame = 30;
ini_time_vector = 10:10*P.n_frame/60:60;

S = zeros(P.n_frame,length(ini_time_vector)-1);
D = zeros(P.n_frame,length(ini_time_vector)-1);
O = zeros(P.n_frame,length(ini_time_vector)-1);
roi_area = zeros(P.n_frame,length(ini_time_vector)-1);

for tt = 1:length(ini_time_vector)-1
    ini_time = find(aux_t == ini_time_vector(tt)*60);
    end_time = find(aux_t == ini_time_vector(tt+1)*60);

    time_lapse = ini_time:10/(P.delta_t):end_time;
    time_lapse = time_lapse(1:end-1);

    sampling = zeros(P.n_frame,P.n*2);
    coeff = zeros(P.n_frame,P.n+1,2);
    coeff2 = zeros(P.n_frame,P.n+1,2);
    coeff_trial = zeros(P.n_frame,P.n+1);
    D_theta = zeros(P.n_frame,P.n*2);
    O_theta = zeros(P.n_frame,P.n*2);
    I_theta = zeros(P.n_frame,P.n*2);

    for l = 1:P.n_frame
%         obtain shape from simulation
        poly = aux_S{time_lapse(l),1};
%         pbtain position of neck center
        Cx = (poly(aux_index3(time_lapse(l),1),1)+poly(aux_index3(time_lapse(l),2),1))/2;
        Cy = poly(aux_index3(time_lapse(l),1),2);
        
%         obtain intersection points
        for jj = 1:2*P.n
            lines = [Cx + P.r*cos(P.theta(jj)), Cy + P.r*sin(P.theta(jj))];
            in = inpolygon(lines(:,1),lines(:,2),poly(:,1),poly(:,2));
            if sum(in)>0
                if sum(in) > 1
                    aux = sqrt((lines(:,1)-Cx).^2 + (lines(:,2)-Cy).^2);
                    aux = aux.*in;
                    sampling(l,jj) = aux(aux == max(aux));
                else
                    sampling(l,jj) = aux(in);
                end
            end
        end
%         obtain Fourier coefficients
        coeff(l,1,1) = sum(sampling(l,:))/(P.n*2);
        for k=1:P.n-1
            coeff(l,k+1,1) =  sum(sampling(l,:).*cos(k*P.theta))./P.n;
            coeff(l,k+1,2) =  sum(sampling(l,:).*sin(k*P.theta))./P.n;
        end
        coeff(l,P.n+1,1) = sum(sampling(l,:).*cos(P.n*P.theta))./(2*P.n);

%         significance test
%         check coefficients with a standard regression coefficient smaller than 0.1 
        a_k = zeros(P.n+1,1);
        b_k = zeros(P.n+1,1);
        s_y = std(sampling(l,:));
        a_k(1) =  coeff(l,1,1)/s_y;
        for k = 1:P.n-1
            a_k(k+1) = coeff(l,k+1,1).*std(cos(k*P.theta))/s_y;
            b_k(k+1) = coeff(l,k+1,2).*std(sin(k*P.theta))/s_y;
        end
        a_k(P.n+1) = coeff(l,P.n+1,1).*std(cos(P.n*P.theta))/s_y;
        aux_a = find(a_k<P.b_min);
        aux_null_a = setdiff(1:(P.n+1),aux_a);
        aux_b = find(b_k<P.b_min);
        aux_null_b = setdiff(1:(P.n+1),aux_b);
        
%         discard the coefficients if they do not significantly alter
%         dROI(theta), with a F-test with level of significance of 0.1
        aux_coeff = zeros(P.n+1,2);
        aux_coeff_null = zeros(P.n+1,2);
        aux_coeff_null(aux_null_a,1) = coeff(l,aux_null_a,1);
        aux_coeff_null(aux_null_b,2) = coeff(l,aux_null_b,2);

        p1 = length(aux_null_a);
        p2 = p1+1;
        for ll=1:length(aux_a)
            aux_coeff(:,:) = aux_coeff_null;
            aux_coeff(aux_a(ll),1) = coeff(l,aux_a(ll),1);
            RSS2 = sum((sampling(l,:)'-IR(1,aux_coeff(:,1)',aux_coeff(:,2)',P)).^2);
            RSS1 = sum((sampling(l,:)'-IR(1,aux_coeff_null(:,1)',aux_coeff_null(:,2)',P)).^2);
            F = ((RSS1-RSS2)/(p2-p1))/(RSS2/(2*P.n-p2));
            F_aux = fpdf(F,p2-p1,2*P.n-p2);
            if F_aux < P.alpha
                p1 = p1+1;
                p2 = p2+1;
                aux_coeff_null = aux_coeff;
            end
        end

        p1 = length(aux_null_b);
        p2 = p1+1;
        for ll=1:length(aux_b)
            aux_coeff(:,:) = aux_coeff_null;
            aux_coeff(aux_b(ll),2) = coeff(l,aux_b(ll),2);
            RSS2 = sum((sampling(l,:)'-IR(1,aux_coeff(:,1)',aux_coeff(:,2)',P)).^2);
            RSS1 = sum((sampling(l,:)'-IR(1,aux_coeff_null(:,1)',aux_coeff_null(:,2)',P)).^2);
            F = ((RSS1-RSS2)/(p2-p1))/(RSS2/(2*P.n-p2));
            F_aux = fpdf(F,p2-p1,2*P.n-p2);
            if F_aux < P.alpha
                p1 = p1+1;
                p2 = p2+1;
                aux_coeff_null = aux_coeff;
            end
        end
        coeff2(l,:,:) = aux_coeff_null;
        coeff_trial(l,:) = ((sqrt(aux_coeff_null(:,1).^2 + aux_coeff_null(:,2).^2))>0);
        [I_theta(l,:),S(l,tt),D_theta(l,:),O_theta(l,:)] = IR(2,coeff2(l,:,1),coeff2(l,:,2),P);
%         Integrate and average D and O
        D(l,tt)=sum(abs(D_theta(l,:))*P.delta_theta)*(100/(S(l,tt)*2*pi));
        O(l,tt)=sum(abs(O_theta(l,1:P.pi_theta)-O_theta(l,2:P.pi_theta+1)))*(100/(S(l,tt)*pi));
   
%         Calculate the are from the sampled ROI
         poly_p = [poly(2:end,:); poly(1,:)];
         roi_area(l,tt) = (1/2).*sum(poly(:,1).*poly_p(:,2)-poly(:,2).*poly_p(:,1));
    end
end
