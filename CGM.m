function [xk,dk,alk,iWk,betak,Hk,tauk] = CGM(x,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,icg,irc,nu)
    tauk = [];
    Hk = [];
    
    n = length(x);
    xk = [x];
    betak = [];
    alk = [];
    iWk = [];
    d = -g(x);
    dk = [d];
    k = 0;
    while norm(g(x)) > epsG && k < kmax  
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
        x = x + al * d;
        xk = [xk,x];
        
        k = k+1;
        
        if icg == 1 %FR
            b = (g(x)'*g(x))/norm(g(xklast))^2;
           betak = [betak,b]; 

        elseif icg == 2 %PR
            b = (g(x)'*(g(x)-g(xklast)))/norm(g(xklast))^2;
            b = max(0,b);
            betak = [betak,b]; 

        end
        
        %restart
        if irc == 0 %RC0
           d = - g(x) + b * d;
        elseif irc == 1 %RC1
           if mod(k,n) == 0
             d = - g(x);
           else
             d = - g(x) + b * d;
           end
         elseif irc == 2 %RC2
             v = abs(g(x)'*g(xklast))/norm(g(x))^2;
             if v >= nu
                 d = - g(x);
             else
                 d = - g(x) + b * d;
             end
        end
       dk = [dk,d];
       
       tauk = [tauk,0];
       Hk = [Hk,0];
    end
end
