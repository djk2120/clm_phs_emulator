function [vwp,qr] = getvwp( q,smp,params )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    k_soil_root = params{1};
    grav1       = params{2};
    lai         = params{3};
    z           = params{4};
    kmax        = params{5};
    p50         = params{6};
    ck          = params{7};
    grav2       = z*1000;
    if sum(k_soil_root)>0
        rootwp = -(q - sum(k_soil_root.*(smp(:)-grav1)))/sum(k_soil_root);
        qr     = k_soil_root.*(smp-rootwp-grav1);
        fx     = plc(rootwp,p50,ck);
        leafwp = rootwp-grav2-q/(lai*fx*kmax/z);
    else
        qr = 0*dz;
    end


    vwp = [rootwp;leafwp];
end

