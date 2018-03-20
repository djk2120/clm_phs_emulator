function f = plc( psi,p50,ck )
    f = 2.^-((psi./p50).^ck);
end