function p=p_normlz(p)

% function normalized_p=p_normlz(p)
% INPUT:
%	p point to be normalized
% OUTPUT:
%	normalized_p normalized point

% (c) DDE-BIFTOOL v. 2.00, 21/10/2001

switch p.kind
  case {'hopf','fold'},
    p.v=p.v/norm(p.v);
  case 'hcli',
    if ~isempty(p.v),
     for k=1:size(p.v,2) 
       p.v(:,k)=p.v(:,k)/norm(p.v(:,k));
     end;
    end;
    if ~isempty(p.w),
      for k=1:size(p.w,2) 
        p.w(:,k)=p.w(:,k)/norm(p.w(:,k));
      end;
    end;
    if ~isempty(p.alpha)
      p.alpha=p.alpha/norm(p.alpha);
    end;
end;

return;
