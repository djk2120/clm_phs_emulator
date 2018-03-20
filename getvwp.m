function [vwp,qr] = getvwp( q,smp,params )
%getvwp Solve for vwp and qr, given a certain q
%   Forced by a given transpiration flux and soil profile, solve for root
%   and leaf water potential. Leaf water potential is used to prognose
%   stress.
%
%   q   , [1]    , transpiration flux (mm/s)
%   smp , [ns,1] , soil potential by soil layer    (mm)
%   vwp , [2,1]  , root and leaf water potential   (mm)
%   qr  , [ns,1] , root water uptake by soil layer (mm/s)

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


    vwp = [leafwp,rootwp];
end

