function f_mem = force_membrane_2D(S,P)
%  S is the vector containing the location of the vertices
% instead of for loops I used vector and matrices, the notation is related
% with the theory S_x = s_i^x, S_x_p = s_{i+1}^x
% pomega_ps = d\Omega / ds
    S_x = S(:,1);
    S_x_p = [S(2:end,1); S(1,1)];
    S_x_p2 = [S(3:end,1); S(1:2,1)];
    S_x_m = [S(end,1); S(1:end-1,1)];
    S_x_m2 = [S(end-1:end,1); S(1:end-2,1)];
    S_y = S(:,2);
    S_y_p = [S(2:end,2);S(1,2)];
    S_y_p2 = [S(3:end,2); S(1:2,2)];
    S_y_m = [S(end,2);S(1:end-1,2)];
    S_y_m2 = [S(end-1:end,2); S(1:end-2,2)];
    
    v = sqrt((S_x-S_x_m).^2 + (S_y-S_y_m).^2);
    v_p = sqrt((S_x_p-S_x).^2 + (S_y_p-S_y).^2);
    v_p2 = sqrt((S_x_p2-S_x_p).^2 + (S_y_p2-S_y_p).^2);
    v_m = sqrt((S_x_m-S_x_m2).^2 + (S_y_m-S_y_m2).^2);
    
    s = (v+v_p)./2;
    s_m = (v_m+v)./2;
    s_p = (v_p + v_p2)./2;
    
    g = ((S_x_p-S_x)./v_p - (S_x-S_x_m)./v).^2 + ...
        ((S_y_p-S_y)./v_p - (S_y-S_y_m)./v).^2;
    g_p = ((S_x_p2-S_x_p)./v_p2 - (S_x_p-S_x)./v_p).^2 + ...
        ((S_y_p2-S_y_p)./v_p2 - (S_y_p-S_y)./v_p).^2;
    g_m = ((S_x-S_x_m)./v - (S_x_m-S_x_m2)./v_m).^2 + ...
        ((S_y-S_y_m)./v - (S_y_m-S_y_m2)./v_m).^2;
    
    pomega_ps = (1/2).*([S_y_p-S_y_m,S_x_m-S_x_p]);
    pSmem_ps = [(S_x-S_x_p)./v_p + (S_x-S_x_m)./v, ...
        (S_y -S_y_p)./v_p + (S_y -S_y_m)./v];
    
   pgk_psk_x = 2*((S_x_p- S_x)./v_p - (S_x-S_x_m)./v).*...
       (((S_x-S_x_p).^2)./(v_p.^3) + ((S_x - S_x_m).^2)./(v.^3) -1./v - 1./v_p)+...
       2*((S_y_p-S_y)./v_p -(S_y-S_y_m)./v).*...
       (((S_y-S_y_p).*(S_x - S_x_p))./(v_p.^3) + ((S_y-S_y_m).*(S_x-S_x_m))./(v.^3));
   pgk_psk_y = 2*((S_y_p- S_y)./v_p - (S_y-S_y_m)./v).*...
       (((S_y-S_y_p).^2)./(v_p.^3) + ((S_y - S_y_m).^2)./(v.^3) -1./v - 1./v_p)+...
       2*((S_x_p-S_x)./v_p -(S_x-S_x_m)./v).*...
       (((S_x-S_x_p).*(S_y - S_y_p))./(v_p.^3) + ((S_x-S_x_m).*(S_y-S_y_m))./(v.^3));
   
   pgkp_psk_x = 2*((S_x_p2-S_x_p)./v_p2 - (S_x_p-S_x)./(v_p)).*...
       (1./v_p - ((S_x_p - S_x).^2)./(v_p.^3)) +...
       2*((S_y_p2-S_y_p)./v_p2 -(S_y_p-S_y)./v_p).*...
       (((S_y_p-S_y).*(S_x-S_x_p))./(v_p.^3)); 
   pgkp_psk_y = 2*((S_y_p2-S_y_p)./v_p2 - (S_y_p-S_y)./(v_p)).*...
       (1./v_p - ((S_y_p - S_y).^2)./(v_p.^3)) +...
       2*((S_x_p2-S_x_p)./v_p2 -(S_x_p-S_x)./v_p).*...
       (((S_x_p-S_x).*(S_y-S_y_p))./(v_p.^3));
   
   pgkm_psk_x = 2*((S_x-S_x_m)./v - (S_x_m-S_x_m2)./v_m).*...
       (1./v - ((S_x-S_x_m).^2)./(v.^3))+...
       2*((S_y - S_y_m)./v - (S_y_m-S_y_m2)./(v_m)).*...
       ((S_y_m-S_y).*(S_x-S_x_m)./(v.^3));
   pgkm_psk_y = 2*((S_y-S_y_m)./v - (S_y_m-S_y_m2)./v_m).*...
       (1./v - ((S_y-S_y_m).^2)./(v.^3))+...
       2*((S_x - S_x_m)./v - (S_x_m-S_x_m2)./(v_m)).*...
       ((S_x_m-S_x).*(S_y-S_y_m)./(v.^3));
   
   pH_ps_x = pgkm_psk_x./s_m - (g_m./(s_m.^2)).*(1/2).*((S_x - S_x_m)./v)+...
       pgk_psk_x./s - (g./(s.^2)).*(1/2).*((S_x-S_x_p)./v_p + (S_x - S_x_m)./v)+...
       pgkp_psk_x./s_p - (g_p./(s_p.^2)).*(1/2).*((S_x - S_x_p)./v_p);
   
   pH_ps_y = pgkm_psk_y./s_m - (g_m./(s_m.^2)).*(1/2).*((S_y - S_y_m)./v)+...
       pgk_psk_y./s - (g./(s.^2)).*(1/2).*((S_y-S_y_p)./v_p + (S_y - S_y_m)./v)+...
       pgkp_psk_y./s_p - (g_p./(s_p.^2)).*(1/2).*((S_y - S_y_p)./v_p);
   
   
   f_mem = -P.P*pomega_ps - P.tau*pSmem_ps - 2*P.kappa*[pH_ps_x, pH_ps_y];

end

