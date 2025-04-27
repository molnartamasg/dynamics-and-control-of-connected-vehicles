function y=LangKobayashi(E,Etau,n,alpha,p,eta,phase,epsilon)
%% r.h.s of Lang-Kobayashi equation
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
Edot=(1+1i*alpha)*n.*E+eta*exp(1i*phase)*Etau;
ndot=epsilon*(p-n-conj(E).*E.*(2*n+1));
y=cat(1,real(Edot),imag(Edot),real(ndot));
end
