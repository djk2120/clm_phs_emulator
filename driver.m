clear
close all

kmax  = 7e-9;
krmax = 1e-9;
lai   = 4;

ck    = 3.95;
p50   = -1e5;


p       = 1;
patchwt = [1];

z     = 25;
soillayers = (0:0.1:2)';
dz    = soillayers(2:end)-soillayers(1:end-1);
soilz = 0.5*(soillayers(1:end-1)+soillayers(2:end));

grav1 = soilz*1000;
ns    = length(soilz);


smp = zeros(ns,1)-10000;

for i=1:1000

fs  = 2.^-((smp/p50).^ck);
k_soil_root = fs.*krmax.*ones(ns,1)./soilz;



q      = 1e-4;
rootwp = -(q - sum(k_soil_root.*(smp(:)-grav1)))/sum(k_soil_root);
qr     = k_soil_root.*(smp-rootwp-grav1);
leafwp = rootwp-q/(lai*kmax/z);


if 3>10*rand
    pp = 3e-4;
else
    pp = 0;
end
smp = bucket( smp,qr,dz,pp);

if mod(i,10)==0
    subplot(1,2,1)
    barh(-1:-1:-ns,smp)
    xlim([-50000,0])
    subplot(1,2,2)
    barh(-1:-1:-ns,qr)
    xlim([-1e-4,3e-4])
    title(i)
    pause(0.1)
end

end

