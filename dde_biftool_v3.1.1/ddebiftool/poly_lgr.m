function p=poly_lgr(t,c)
%% values of lagrange polynomials through t at c
% function p=poly_lgr(t,c);
% INPUT:
%	t lagrange points in R^m+1
% 	c evaluation point(s) in R 
% OUTPUT:
%	p values of lagrange polynomials through t at c
%

% (c) DDE-BIFTOOL v. 3.1.1(19), 11/04/2014
%
% 
%
%%
m=length(t)-1;
nc=length(c);

% compute p:
p=ones(m+1,nc);
for j=1:m+1
    for k=1:m+1
        if k~=j
            p(j,:)=p(j,:).*(c-t(k))/(t(j)-t(k));
        end
    end
end
% reshape to keep compatibility with original version
if nc==1
    p=p';
end
end
