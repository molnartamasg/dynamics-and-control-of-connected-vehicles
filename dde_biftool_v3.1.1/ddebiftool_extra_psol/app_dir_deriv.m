function dy=app_dir_deriv(f,x0,deviation,ind,h)
%% 2nd order approx directional derivative in ith column of x0
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
xd=x0;
xd(:,ind,:)=x0(:,ind,:)+h*deviation;
dy1=f(xd);
xd(:,ind,:)=x0(:,ind,:)-h*deviation;
dy2=f(xd);
dy=(dy1-dy2)/(2*h);
end