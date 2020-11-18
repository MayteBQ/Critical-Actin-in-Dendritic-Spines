function S_int = solve_system_RK(S,aux_fil,actin_dyn,delta_t,P)
% Runge Kutta 4 method
    if actin_dyn == 0 
        k1 = delta_t*P.zeta*(force_membrane_2D(S,P));
        k2 = delta_t*P.zeta*(force_membrane_2D(S+k1./2,P));
        k3 = delta_t*P.zeta*(force_membrane_2D(S+k2./2,P));
        k4 = delta_t*P.zeta*(force_membrane_2D(S+k3,P));
    else
        k1 = delta_t*P.zeta*(force_membrane_2D(S,P)-f_fil_2D(S,aux_fil,P));
        k2 = delta_t*P.zeta*(force_membrane_2D(S+k1./2,P)-f_fil_2D(S+k1./2,aux_fil,P));
        k3 = delta_t*P.zeta*(force_membrane_2D(S+k2./2,P)-f_fil_2D(S+k2./2,aux_fil,P));
        k4 = delta_t*P.zeta*(force_membrane_2D(S+k3,P)-f_fil_2D(S+k3,aux_fil,P));
    end
    S_int = S + (k1+2*k2+2*k3+k4)/6;
end