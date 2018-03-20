function [vwp,q,qr] = phs(x)
%phs Solves for vwp and qr, subject to parameters, smp, and qmax

kmax  = x{1};
krmax = x{2};
lai   = x{3};
ck    = x{4};
p50   = x{5};
z     = x{6};        %tree height
soillayers = x{7};   %soil layer boundaries         [ns+1,1] (m)
smp   = x{8};        %soil water potential          [ns,1]   (mm)
qmax  = x{9};        %maximum transpiration flux    [1]      (mm/s)
rai   = x{10};       %root area index by soil layer [ns,1]
hk    = x{11};       %soil hydraulic conductivity   [ns,1]   (mm/s)
soilz = 0.5*(soillayers(1:end-1)+soillayers(2:end));
grav1 = soilz*1000;  %soil layer gravitational potential to surface

% solve for k_soil_root
fs =plc(smp,p50,ck);
root_cond = fs.*rai*krmax./(soilz+0.25);
r_soil = max(rai)./rai; %this is just a rough placeholder
soil_cond = hk./r_soil;
k_soil_root = (1./root_cond+1./soil_cond).^-1;

% package parameters
params = cell(7,1);
params{1}=k_soil_root;
params{2}=grav1;
params{3}=lai;
params{4}=z ;
params{5}=kmax;
params{6}=p50;
params{7}=ck;

%iterate to match transpiration flux (q2)
% with root water uptake (q1)
% starting at qmax and stepping towards 0
% note this is not quite how it's solved in CLM
% CLM uses Newton's method
% but this is easier to code up and should work for a first pass
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
        q1 = sum(qr);                 %root water uptake
        q2 = qmax*plc(vwp(1),p50,ck); %transpiration

        if abs(q1-q2)<0.0001*q1
            go = 0;
        elseif q1<1e-10
            go = 0;
        elseif q2>q1
            q = q+step;
            step = step/2;
        end
        
        if ct>50
            go = 0;
            'Not converging well'
        end
    end
    
end
    q = q1;
    
end