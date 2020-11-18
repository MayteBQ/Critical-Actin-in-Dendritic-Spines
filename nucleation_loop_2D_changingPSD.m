function [ff_length,ini_fil,a_points,line,filaments] = nucleation_loop_2D_changingPSD(t,ff_length,filaments,S,P)
    ini_fil = P.ini_fil;
    a_points = P.a_points;
    line = P.line;
 
    if rand < P.delta_t*P.tau_n
        [ff_length,ini_fil,a_points,line,aux_fil] = add_B_2D_loop_changingPSD(t,ff_length,S,P);
         filaments = [filaments;aux_fil];
    end
end