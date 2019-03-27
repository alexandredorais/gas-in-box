function[P_exp, P_th]=calculePression(deltap_axe,S,deltat, N, T, V)
    % pression gaz parfaits
    kB = 1.3806e-23;
    P_th = N*kB*T/V;
    % pression experimentale
    VecP=deltap_axe./(deltat.*S);
    P_exp=sum(VecP);
end