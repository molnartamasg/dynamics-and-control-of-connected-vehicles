function y=VAopX(A,x,op)
%% vectorized matrix multiplication/linear system solving
%
% function y=VAopX(A,x,op)
%
%applies matrix operation (y=A*x or y=A\x) to higher dim arrays A(1:na1,1:na2,k),
%and x(1:nx1,1:nx2,k) or x(1:nx1,k), returning 
%y(:,:,k)=A(:,:,k) op x(:,:,k) or
%y(:,k)=A(:,:,k) op x(:,k) or
%where k can be a higher dimensional index, all dimensions of A and x have
%to match: no "scalar" expansion is done except for the second dimension of x.
%
%op can be a symbol: op='*', '\'.
%
%Relies on sparse_blkdiag, works also for larger numbers of
%dimensions. May be more efficient than loops if na1,na2 are small but na3 is
%large.
%
%
% (c) DDE-BIFTOOL v. 3.1.1(19), 11/04/2014
%
%%
dima=size(A);
dimx=size(x);
ndimx=length(dimx);
ndima=length(dima);
isvec=0;
if ndimx==ndima-1
    %2nd Arg is vector
    isvec=1;
    x=reshape(x,[dimx(1),1,dimx(2:end)]);
    dimx=[dimx(1),1,dimx(2:end)];
    ndimx=ndimx+1;
end
nveca=prod(dima(3:end));
nvecx=prod(dimx(3:end));
Ar=reshape(A,[dima(1),dima(2),nveca]);
xr=reshape(x,[dimx(1),dimx(2),nvecx]);
if op=='*'
    xr=permute(xr,[2,1,3]);
    xr=reshape(xr,[dimx(2),dimx(1)*nvecx]);
    xr=xr';
    yr=sparse_blkdiag(Ar)*xr;
    yr=reshape(yr',[dimx(2),dima(1),dimx(3:end)]);
    y=permute(yr,[2,1,3:ndimx]);
elseif op=='\'
    xr=permute(xr,[2,1,3]);
    xr=reshape(xr,[dimx(2),dimx(1)*nvecx]);
    xr=xr';
    yr=sparse_blkdiag(Ar)\xr;
    yr=reshape(yr',[dimx(2),dimx(1),dimx(3:end)]);
    y=permute(yr,[2,1,3:ndimx]);
else
    error('VAopX: operation "%s" not implemented\n',op);
end
if isvec
    y=reshape(y,[dimx(1),dimx(3:end)]);
end
end
