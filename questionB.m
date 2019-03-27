clear all; clc

N=100; % nbr de particules
uma=1.66e-27;
mH2 = 2*uma; % masse hydrogene
rH2 = 0.120e-9; % rayon hydrogene
global mN2
mN2 = 2*14*uma;
rN2 = 0.155e-9;
L = 15e-9; % taille boite
S = 6*L^2; % surface boite
V = L^3; % volume boite
T = 300;
repetitions = 2000; % nombre de pas de temps consideres 

particules = zeros(N+1,7); % on initialise la matrice de proprietes des particules
% les colonnes de 'particules' sont en ordre:
% 1 - masse
% 2 - position en x
% 3 - position en y
% 4 - position en z
% 5 - vitesse en x
% 6 - vitesse en y
% 7 - vitesse en z
particules(1:end-1,1) = mH2;
particules(1:end-1,2) = L*rand(N,1) - L/2;
particules(1:end-1,3) = L*rand(N,1) - L/2;
particules(1:end-1,4) = L*rand(N,1) - L/2;
v_rand = maxBoltz(mH2,T,N); % genere selon maxwell-boltzman N normes de v
theta_rand = pi*rand(N,1);
phi_rand = 2*pi*rand(N,1);
vx = v_rand.*sin(theta_rand).*cos(phi_rand);
vy = v_rand.*sin(theta_rand).*sin(phi_rand);
vz = v_rand.*cos(theta_rand);
particules(1:end-1,5) = vx;
particules(1:end-1,6) = vy;
particules(1:end-1,7) = vz;
particules(end,1) = mN2;

% si toutes les particules donnaient leur energie a la meme particule, et
% que cette derniere etait la seule avec une vitesse, quelle serait sa
% vitesse? : v_th_max
v_calc_max = sqrt(sum(v_rand.^2));
% et si toutes les particules avaient la meme energie, quelle serait leur
% vitesse?
v_exp_moy = v_calc_max/sqrt(N);
% ainsi, un bon increment de temps pour la simulation est delta_t
delta_t = rH2/(2*v_calc_max)*5;
temps_total = delta_t*repetitions;

H=figure;
scatter3(particules(1:end-1,2),particules(1:end-1,3),particules(1:end-1,4),'filled','k')
hold on
scatter3(particules(end,2),particules(end,3),particules(end,4),'filled','r')
hold off
xlim([-L/2 L/2])
ylim([-L/2 L/2])
zlim([-L/2 L/2])
drawnow

% premier pas de temps, on update les positions
% j=2 : update en x; j=3 : update en y; j=4 : update en z
particules(:,2:4) = particules(:,2:4) + particules(:,5:7).*delta_t;
scatter3(particules(1:end-1,2),particules(1:end-1,3),particules(1:end-1,4),'filled','k')
hold on
scatter3(particules(end,2),particules(end,3),particules(end,4),'filled','r')
hold off
xlim([-L/2 L/2])
ylim([-L/2 L/2])
zlim([-L/2 L/2])
drawnow

global nombre_collisions
nombre_collisions=0;
variation_p = [];
i=0;
tic
% tous les autres pas de temps
while nombre_collisions < 10
    i=i+1;
    % collisions avec mur (conditions : 1) pres du mur, 2) vitesse vers le
    % mur)
    condition_murXplus = abs(L/2-particules(:,2)) <= rH2 & particules(:,5) > 0;
    condition_murXmoins = abs(-L/2-particules(:,2)) <= rH2 & particules(:,5) < 0;
    condition_murX = condition_murXplus | condition_murXmoins;
    condition_murYplus = abs(L/2-particules(:,3)) <= rH2 & particules(:,6) > 0;
    condition_murYmoins = abs(-L/2-particules(:,3)) <= rH2 & particules(:,6) < 0;
    condition_murY = condition_murYplus | condition_murYmoins;
    condition_murZplus = abs(L/2-particules(:,4)) <= rH2 & particules(:,7) > 0;
    condition_murZmoins = abs(-L/2-particules(:,4)) <= rH2 & particules(:,7) < 0;
    condition_murZ = condition_murZplus | condition_murZmoins;
    % on sauve les variations de qte de mouvement
    tempX = abs(2.*particules(condition_murX,1).*particules(condition_murX,5)); % 2*m*|v_x|
    tempY = abs(2.*particules(condition_murY,1).*particules(condition_murY,6)); % 2*m*|v_y|
    tempZ = abs(2.*particules(condition_murZ,1).*particules(condition_murZ,7)); % 2*m*|v_z|
    variation_p = [variation_p; tempX; tempY; tempZ];
    clear tempX tempY tempZ
    % on modifie les vitesses
    particules(condition_murX,5) = - particules(condition_murX,5);
    particules(condition_murY,6) = - particules(condition_murY,6);
    particules(condition_murZ,7) = - particules(condition_murZ,7);
    
    % collisions entre particules
    for n = 1:N+1 % premiere particule
        for m = n:N+1 % deuxieme particule
            if (n < N+1) && (m < N+1)
                if norm(particules(n,2:4) - particules(m,2:4)) <= rH2 + rH2
                    [vf1, vf2] = collisionB(particules(n,5:7), particules(m,5:7), particules(n,2:4), particules(m,2:4), particules(n,1), particules(m,1));
                    % dans la fonction 'collision' on verifie des conditions
                    %   sur la vitesse pour s'assurer d'avoir collision, donc 
                    %   on incremente nombre_collisions dans la fonction
                    % on modifie les vitesses
                    particules(n,5:7) = vf1;
                    particules(m,5:7) = vf2;
                end
            elseif m == N+1
                if norm(particules(n,2:4) - particules(m,2:4)) <= rH2 + rN2
                    [vf1, vf2] = collisionB(particules(n,5:7), particules(m,5:7), particules(n,2:4), particules(m,2:4), particules(n,1), particules(m,1));
                    % dans la fonction 'collision' on verifie des conditions
                    %   sur la vitesse pour s'assurer d'avoir collision, donc 
                    %   on incremente nombre_collisions dans la fonction
                    % on modifie les vitesses
                    particules(n,5:7) = vf1;
                    particules(m,5:7) = vf2;
                end
            end
        end
    end
    
    % update les positions des molecules
    particules(:,2:4) = particules(:,2:4) + particules(:,5:7).*delta_t;
    
    % si les particules sont a l'exterieur de la boite par erreur, on les
    % ramene dedans
    exterieur = abs(particules(:,2:4))> L/2; 
    particules(exterieur(:,1),2) = L/2 - abs(particules(exterieur(:,1),2) - L/2);
    particules(exterieur(:,2),3) = L/2 - abs(particules(exterieur(:,2),3) - L/2);
    particules(exterieur(:,3),4) = L/2 - abs(particules(exterieur(:,3),4) - L/2);
    % et on inverse leur vitesse
    particules(exterieur(:,1),5) = -particules(exterieur(:,1),5);
    particules(exterieur(:,2),6) = -particules(exterieur(:,2),6);
    particules(exterieur(:,3),7) = -particules(exterieur(:,3),7);
    % il faut aussi sauver les variations de qte de mouvement, car c'est
    % equivalent a une collision sur le mur
    tempX = abs(2*particules(exterieur(:,1),1).*particules(exterieur(:,1),5)); % 2*m*|v_x|
    tempY = abs(2*particules(exterieur(:,2),1).*particules(exterieur(:,2),6)); % 2*m*|v_y|
    tempZ = abs(2*particules(exterieur(:,3),1).*particules(exterieur(:,3),7)); % 2*m*|v_z|
    variation_p = [variation_p; tempX; tempY; tempZ];
    clear tempX tempY tempZ
    
    
    % dessine
%      scatter3(particules(1:end-1,2),particules(1:end-1,3),particules(1:end-1,4),'filled','k')
%      hold on
%      scatter3(particules(end,2),particules(end,3),particules(end,4),'filled','r')
%      hold off
%      xlim([-L/2 L/2])
%      ylim([-L/2 L/2])
%      zlim([-L/2 L/2])
%      drawnow
end

% déplacement total de la molécule d'azote
Rtot = norm(particules(end,2:4));



