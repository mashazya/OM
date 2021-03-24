function [al,iWout] = bls(xk,dk,f,g,almax,almin,rho,c1,c2,iW)
    phi = @(al) f(xk+al*dk);
    dphi = @(al) g(xk+al*dk)'*dk;
    WC1 = @(al) phi(al) <= phi(0) + c1*dphi(0)*al;
    WC2 = @(al) dphi(al) >= c2*dphi(0);
    WC = @(al) WC1(al) & WC2(al);
    SWC2 = @(al) abs(dphi(al)) <= c2* abs(dphi(0));
    SWC = @(al) WC1(al) & SWC2(al);
    ak = almax;
    iWout = 0;
    while ak >= almin
       if WC(ak) && iW == 1
           iWout = 2;
           break;
       elseif SWC(ak) && iW == 2
           iWout = 3;
           break;
       elseif WC1(ak)
           iWout = 1;
       end
       ak= rho*ak;
    end
    al = ak;
end
