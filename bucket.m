function [psoil1,sm1,ee,hk] = bucket( psoil0,qr,dz,p,evapflag )
%bucket Simplified bucket model
%   inputs
%     psoil0, original soil water potential [mm H2O]
%     q     , ET                            [mm/s]
%     dz    , layer thickness               [m]
%     p     , Precip                        [mm/s]
%     evapflag  , turn soil evap on/off
qr = qr/1e3;  %convert to meters/s
p  = p/1e3;
dt = 1800;   %timestep = 30minutes = 1800seconds


b = 5;
psat  = -300;
smsat = 0.5;
hksat = 1e-3; %mm/s


% calc soil moisture
sm0   = smsat.*(psoil0/psat).^(-1/b);


%check for available water
if sum(sm0.*dz)<sum(qr.*dt)
    error('Error soils too dry')
end

%subtract transpiration
sm1   = sm0-qr*dt./dz;

%add in precip water
if p*dt>(smsat-sm1(1))*dz(1)
    i  = 0;
    go = 1;
    while go
    i = i+1;
    
    if i>length(psoil0)
        go=0;
    else
        space_for_water = (smsat-sm1(i))*dz(i);
        if space_for_water>p*dt
            sm1(i)=sm1(i)+p*dt./dz(1);
            go = 0;
        else
            sm1(i)=smsat;
            p     = p-(p*dt-space_for_water)/dt;
        end
    end
    end
else
    sm1(1)= sm1(1)+p*dt./dz(1);
end

%subtract soil evap
if evapflag
    if p<1e-6
        ee = 0.1*(exp((sm1(1)./0.5).^2)-1);  %mm/halfhour (very crude)
        if ee<0.01
            ee=0;
        end
        sm1(1) = sm1(1) - ee/1000./dz(1);
    end
else
    ee=0;
end



psoil1 = psat.*(sm1/smsat).^-b;
hk     = hksat*(sm1/smsat).^(2*b+3);

end