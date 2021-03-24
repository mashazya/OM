function [xk,dk,alk,iWk,betak,Hk,tauk] = NM(x,g,h,epsG,kmax)
    tauk = [];
    betak = [];
    Hk = [];
    xk = [x];
    alk = [];
    iWk = [];
    dk = [];
    k = 0;
    while norm(g(x)) > epsG && k < kmax-1  
        d = - h(x)^-1 * g(x);
        dk = [dk,d];
        
        x = x + d;
        xk = [xk,x];
        
        k = k+1;
        
        tauk = [tauk,0];
        iWk = [iWk,4];
        alk = [alk,1];
    end
end
