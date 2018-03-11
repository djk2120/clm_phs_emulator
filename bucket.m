function psoil1 = bucket( psoil0,qr,dz,p )
%bucket Simplified bucket model
%   inputs
%     psoil0, original soil water potential [MPa]
%     q     , ET                            [mm/s]
%     dz    , layer thickness               [m]

qr = qr/1e3;  %convert to meters/s
p  = p/1e3;
dt = 1800;   %timestep = 30minutes = 1800seconds


b = 5;
psat  = -300;
smsat = 0.5;

sm0   = smsat.*(psoil0/psat).^(-1/b);
sm1   = sm0-qr*dt./dz;

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



psoil1 = psat.*(sm1/smsat).^-b;

end