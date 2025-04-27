function hcli=p_tohcli(funcs,point)
%% convert point to connecting orbit
% INPUT:
%     funcs problem functions
%     point a periodic solution near a homoclinic solution
%           alternatively an initial point in a hcli structure,
%           where a good starting guess for the profile and steady
%           states are available
% OUTPUT:
%     hcli a starting value to compute the exact homoclinic or
%     heteroclinic solution  

% (c) DDE-BIFTOOL v. 3.1.1(19), 11/04/2014
%
% 
%
%%
sys_tau=funcs.sys_tau;
sys_deri=funcs.sys_deri;

    if mod(length(point.mesh),point.degree)~=1,
      err=[length(point.mesh) point.degree];
      error('P_TOHCLI: psol does not contain L=%d intervals of m=% points!',...
          err(1),err(2));
    end;
    
    hcli.kind='hcli';
    hcli.parameter=point.parameter;
    hcli.mesh=point.mesh;
    hcli.degree=point.degree;
    
    switch point.kind,
     
     case 'psol',
      ntst=size(point.profile,2);   
      test=NaN(1,ntst-1);
      for i=1:ntst-1
          test(i)=norm(point.profile(:,i)-point.profile(:,i+1));
      end;
      [minval, pos]=min(abs(test)); %#ok<ASGLU>
      stst.kind='stst';
      stst.parameter=hcli.parameter;
      stst.x=point.profile(:,pos);
      x_profile=NaN(1,ntst);
      for i=1:size(point.profile,2)
          x_profile(1,i)=norm(point.profile(:,i)-stst.x);
      end;
      
      [peak, peak_pos]=max(x_profile); %#ok<ASGLU>
      [hole, hole_pos]=min(x_profile); %#ok<ASGLU>
      left_part=point.profile(:,1:peak_pos);
      right_part=point.profile(:,peak_pos+1:end);
      hole_begin=hole_pos-mod(hole_pos,point.degree)+1;
      hole_end=hole_begin+point.degree;
      
      if hole_pos<peak_pos,
          right_part=[right_part left_part(:,2:hole_begin)];
          left_part=left_part(:,hole_end:end);
          hcli.mesh=[hcli.mesh(hole_end:end) ...
              (hcli.mesh(2:hole_begin)+1)];
      else
          left_part=[right_part(:,hole_end-peak_pos:end-1) left_part];
          right_part=right_part(:,1:hole_begin-peak_pos);
          hcli.mesh= ...
              [(hcli.mesh(hole_end:end-1)-1)...
              hcli.mesh(1:hole_begin)];
      end
      
      nb_of_points=length(hcli.mesh);
      rest=mod(nb_of_points,point.degree);
      hcli.profile=[left_part right_part];
      
      if rest>1,
          hcli.profile=point.profile(:,1+floor((rest-1)/2):end-ceil((rest-1)/2));
          hcli.mesh=hcli.mesh(1+floor((rest-1)/2):end-ceil((rest-1)/2));
      end
      if rest==0,
          rest=point.degree;
          hcli.profile=point.profile(:,1+floor((rest-1)/2):end-ceil((rest-1)/2));
          hcli.mesh=hcli.mesh(1+floor((rest-1)/2):end-ceil((rest-1)/2));
      end
      hcli.mesh=hcli.mesh-hcli.mesh(1);
      hcli.period=point.period*hcli.mesh(end);
      hcli.mesh=hcli.mesh/hcli.mesh(end);
      hcli.x1=stst.x;
      hcli.x2=stst.x;
      stst1=stst;
      stst2=stst;
     case 'hcli',
      hcli=point;
      stst1.kind='stst';
      stst1.parameter=point.parameter;
      stst1.x=point.x1;
      stst2=stst1;
      stst2.x=point.x2;
     otherwise,
      error(['P_TOHCLI: not a valid conversion for other than psol' ...
	     ' or hcli type points']);
    end;
    
    m=df_mthod(funcs,'stst');
    stst1.stability=p_stabil(funcs,stst1,m.stability);
    
    if isempty(stst1.stability.l1) || max(real(stst1.stability.l1))<0
      error('P_TOHCLI: no unstable eigenmodes found');
    end
    lambda=stst1.stability.l1(:);
    lambda=lambda(real(lambda)>0);
    
    hcli.lambda_v=lambda;
    
    if funcs.tp_del==0
      tau=point.parameter(sys_tau());
      n_tau=length(tau);
    else
      error('P_TOHCLI: computing connected orbits is not implemented for equations with state-dependent delays');
    end;
    
    
    n=length(hcli.profile(:,1));
    if n==1,
      v=ones(1,length(lambda));
    else
      v=NaN(n,length(lambda));  
      for i=1:length(lambda)
          delta=eye(n)*lambda(i);
          xx=stst1.x(:,ones(n_tau+1));
          delta=delta-sys_deri(xx,hcli.parameter,0,[],[]);
          for t=1:n_tau
              delta=delta-sys_deri(xx,hcli.parameter,t,[],[])*exp(-lambda(i)*tau(t));
          end
          [eigvec,eigval]=eig(delta);
          [minval pos]=min(abs(diag(eigval))); %#ok<ASGLU>
          v(:,i)=eigvec(:,pos);
      end
    end
    
   stst2.stability=p_stabil(funcs,stst2,m.stability);
   
   lambda=stst2.stability.l1(:);
   lambda=lambda(real(lambda)>0); 
   
   hcli.lambda_w=lambda;
    
   if n==1
       w=ones(1,length(lambda));
   else
       w=NaN(n,length(lambda));
       for i=1:length(lambda),
           delta=eye(n)*lambda(i);
           xx=stst2.x(:,ones(n_tau+1,1));
           delta=delta-sys_deri(xx,hcli.parameter,0,[],[]);
           for t=1:n_tau
               delta=delta-sys_deri(xx,hcli.parameter,t,[],[])*exp(-lambda(i)*tau(t));
           end
           delta=delta';
           [eigvec,eigval]=eig(delta);
           [minval pos]=min(abs(diag(eigval))); %#ok<ASGLU>
           w(:,i)=eigvec(:,pos);
       end
   end
    
   hcli.v=v;
   hcli.w=w;
    
   hcli.alpha=hcli.v\(hcli.profile(:,1)-hcli.x1);
   hcli.epsilon=norm(hcli.alpha);
   hcli.alpha=hcli.alpha/hcli.epsilon;
end
    
