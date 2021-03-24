function [xk,dk,alk,iWk,betak,Hk,tauk] = MNM(x,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW)
    tauk = [];
    betak = [];
    Hk = [];
    xk = [x];
    alk = [];
    iWk = [];
    dk = [];
    k = 0;
    while norm(g(x)) > epsG && k < kmax-1
        d = -h(x)^-1 * g(x);
        dk = [dk,d];
        %find alk
        if iW == 0 %ELS
            %g(x) = Qx - b
            Q = h(x);
            al = -(g(x)'*d)/(d'*Q*d); %direccio del metode o d = -g(x) ??
            iWk = [iWk, -1];
        elseif iW == 1 || iW == 2 %BLS WC or SWC
            [al,iWx] = bls(x,d,f,g,almax,almin,rho,c1,c2,iW);
            iWk = [iWk, iWx];
            alk = [alk,al];
        end
        
        x = x+al*d;
        xk = [xk,x];
        
        k = k+1;
        
        tauk = [tauk,0];
    end
end
