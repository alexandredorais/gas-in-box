function[lambdaexp, lambdatheo]=libreParcoursMoy(P,v_moy,moydeltat)
    % libre parcours moyen experimental
    lambdaexp=v_moy.*moydeltat;
    % libre parcours moyen theorique
    kB = 1.3806e-23;
    T=300;
    sigmaeff=0.27*(1e-9)^2;
    lambdatheo=kB.*T./(P.*sigmaeff);
end