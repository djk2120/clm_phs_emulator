function [vwp,q,qr] = phs(x)

kmax  = x{1};
krmax = x{2};
lai   = x{3};
ck    = x{4};
p50   = x{5};
z     = x{6};
soillayers = x{7};
smp   = x{8};
qmax  = x{9};
rai   = x{10};
hk    = x{11};
soilz = 0.5*(soillayers(1:end-1)+soillayers(2:end));
grav1 = soilz*1000;

fs =plc(smp,p50,ck);
root_cond = fs.*rai*krmax./(soilz+0.25);

r_soil = max(rai)./rai; %this is not robust
soil_cond = hk./r_soil;
k_soil_root = (1./root_cond+1./soil_cond).^-1;

params = cell(7,1);
params{1}=k_soil_root;
params{2}=grav1;
params{3}=lai;
params{4}=z ;
params{5}=kmax;
params{6}=p50;
params{7}=ck;


if qmax ==0
    [vwp,qr] = getvwp( qmax,smp,params );
    q1 = sum(qr);
else
    go =1;
    ct = 0;
    q  = qmax;
    step = qmax/10.5;
    while go
        ct = ct+1;
        if step>q
            step = q/5;
        end
        q = q-step;
        
        [vwp,qr] = getvwp( q,smp,params );
        
        q1 = sum(qr);
        q2 = qmax*plc(vwp(2),p50,ck);
        

        
        if q2>q1
            q = q+step;
            step = step/2;
        end
        
        
        if ct>50
            go = 0;
        end
    end
    
end
    q = q1;
    
end