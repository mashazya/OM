function [xk,dk,alk,iWk,betak,Hk,tauk] = BFGS(x,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW)
    tauk = [];
    
    n = length(x);
    I = eye(n);
    H = I;
    Hk(:,:,1) = H;
    xk = [x];
    betak = [];
    alk = [];
    iWk = [];
    dk = [];
    k = 0;
    while norm(g(x)) > epsG && k < kmax
        d = -H*g(x);
        dk = [dk,d];
        %find alk
        if iW == 0 %ELS
            %g(x) = Qx - b
            Q = h(x);
            al = -(g(x)'*d)/(d'*Q*d);  %direccio del metode o d = -g(x) ??
            iWk = [iWk, -1];
        elseif iW == 1 || iW == 2 %BLS WC or SWC
            [al,iWx] = bls(x,d,f,g,almax,almin,rho,c1,c2,iW);
            iWk = [iWk, iWx];
        end
        alk = [alk,al];
        
        xklast = x;
        x = x + al*d;
        xk = [xk,x];
        
        sk = x - xklast;
        yk = g(x) - g(xklast);
        pk = 1/(yk'*sk);
        
        H = (I-(pk*sk*yk'))*H*(I-(pk*yk*sk'))+pk*sk*sk';
        k = k+1;
        Hk(:,:,k+1) = H;
       
    end
end 