function f = plc( psi,p50,ck )
%plc Returns f, fraction of maximum conductance
    f = 2.^-((psi./p50).^ck);
end