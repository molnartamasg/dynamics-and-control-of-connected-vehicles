function dxdt = OVMNoDelayRing(~,x,alpha,beta,hst,hgo,vmax)
    % headway and followers' velocity
    h=x(1:2:end);       
    vF=x(2:2:end);
    % leaders' velocity (each follower is a leader for the subsequent one)
    vL=[vF(end);...
        vF(1:end-1)];	
    % range policy
    Fh=vmax.*(h-hst)./(hgo-hst);
    Vh=0+(hst<h & h<hgo).*Fh+(hgo<=h).*vmax;
    % optimal velocity model without driver reaction time
    dxdt = [(vL-vF).';(alpha.*(Vh-vF)+beta.*(vL-vF)).'];
    dxdt = dxdt(:);
end
