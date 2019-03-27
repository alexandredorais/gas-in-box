function [v] = maxBoltz(m,T,nb)
    kB = 1.3806e-23;
    vmax = sqrt(2*kB*T/m);
    fmax = 4*exp(-1)*sqrt(m/(2*pi*kB*T));
    vi = linspace(0,10*vmax,100000);
    fi = (m/(2*pi*kB*T))^(3/2)*4*pi*vi.^2.*exp(-m.*vi.^2/(2*kB*T));
    occurences = round(200*fi/fmax);
    vpick = [];
    for i = 1:length(vi)
        temp(1:occurences(i),1) = vi(i);
        vpick = [vpick; temp];
        clear temp
    end
    maxN = length(vpick);
    choisis = randi([1,maxN],1,nb);
    v = vpick(choisis);
end