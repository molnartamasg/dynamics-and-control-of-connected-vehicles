%% Range policy functions
function [V,dV,ddV]=RP(RPtype,vmax,hgo,hst)
    % linear range policy
    if strcmp(RPtype,'lin')
        V=@(h)vmax*(hgo<=h)...
            + vmax/(hgo-hst)*(h-hst).*(hst<h & h<hgo);
        dV=@(h)vmax/(hgo-hst).*(hst<h & h<hgo);
        ddV=@(h)0*h;
        
    % quadratic range policy
    elseif strcmp(RPtype,'quad')
        V=@(h)vmax*(hgo<=h)...
            + vmax*(2*hgo-hst-h).*(h-hst)/(hgo-hst)^2.*(hst<h & h<hgo);
        dV=@(h)vmax*2*(hgo-h)/(hgo-hst)^2.*(hst<h & h<hgo);
        ddV=@(h)-vmax*2/(hgo-hst)^2.*(hst<h & h<hgo);
        
    % cubic range policy
    elseif strcmp(RPtype,'cube')
        V=@(h)vmax*(hgo<=h)...
            + vmax*(3*hgo-hst-2*h).*(h-hst).^2/(hgo-hst)^3.*(hst<h & h<hgo);
        dV=@(h)vmax*6*(hgo-h).*(h-hst)/(hgo-hst)^3.*(hst<h & h<hgo);
        ddV=@(h)vmax*6*(hgo+hst-2*h)/(hgo-hst)^3.*(hst<h & h<hgo);
        
    % cosine range policy
    elseif strcmp(RPtype,'cos')
        V=@(h)vmax*(hgo<=h)...
            + vmax/2*(1-cos(pi*(h-hst)/(hgo-hst))).*(hst<h & h<hgo);
        dV=@(h)vmax/2*sin(pi*(h-hst)/(hgo-hst))*pi/(hgo-hst).*(hst<h & h<hgo);
        ddV=@(h)vmax/2*cos(pi*(h-hst)/(hgo-hst))*pi^2/(hgo-hst)^2.*(hst<h & h<hgo);
    end
end