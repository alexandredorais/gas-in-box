function [v1f, v2f] = collisionB(v1,v2,r1,r2,m1,m2)
global nombre_collisions
global mN2
    %r1,r2,v1,v2 sont des vecteurs a 3 composantes
    u = r2-r1;
    v_para1 = dot(v1,u)/dot(u,u)*u;
    v_para1norm = dot(v1,u)/norm(u); % positif ou negatif
    v_perp1 = v1-v_para1;
    v_para2 = dot(v2,u)/dot(u,u)*u;
    v_para2norm = dot(v2,u)/norm(u); % positif ou negatif
    v_perp2 = v2-v_para2;
    
    condition1 = (dot(v_para1,u) >= 0 & dot(v_para2,u) <= 0);
    condition2 = (dot(v_para1,u) <= 0 & dot(v_para2,u) <= 0 & norm(v_para2) > norm(v_para1) );
    condition3 = (dot(v_para1,u) >= 0 & dot(v_para2,u) >= 0 & norm(v_para1) > norm(v_para2) );
    
    if any([condition1, condition2, condition3]) % collision
        if (m1 == mN2) || (m2 == mN2)
            nombre_collisions = nombre_collisions + 1;
            toc
        end
        matrice = [m1,m2; 1,-1];
        vecteur_i = [m1*v_para1norm + m2*v_para2norm; v_para2norm - v_para1norm];
        v_f = matrice\vecteur_i;
        vf_para1 = v_f(1)*u/norm(u);
        vf_para2 = v_f(2)*u/norm(u);
        v1f = vf_para1 + v_perp1;
        v2f = vf_para2 + v_perp2;
    else  % pas de collision
        v1f = v1;
        v2f = v2;
    end
end