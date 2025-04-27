function hopf=p_tohopf(funcs,point,freqs)
%% convert point to Hopf bifurcation point
% function hopf_point=p_tohopf(funcs,point {,freqs})
% INPUT:
%   funcs problem functions
%	point with stability information 
%   optional freqs: frequency to be excluded from consideration (or
%   used otherwise, depending on input point type)
% OUTPUT:
%	hopf_point uncorrected starting guess for hopf point
%
% If point is of type 'hoho' freqs can be string 'omega1' or 'omega2': then
% omega1 or omega2 will be selected as Hopf frequency, if freqs is value,
% it will be excluded
%
%
% (c) DDE-BIFTOOL v. 3.1.1(80), 02/01/2015
%
%%
if nargin<3
    freqs=[];
end

hopf.kind='hopf';
hopf.parameter=point.parameter;

switch point.kind
    case {'stst','fold','hopf'}
        if ~isfield(point,'stability') || isempty(point.stability)
            error('P_TOHOPF: point does not contain stability information!');
        end
        l1=point.stability.l1;
        if isempty(l1)
            error('P_TOHOPF: point does not contain stability information l1!');
        end
        if strcmp(point.kind,'hopf')
            % remove known imaginary pair
            freqs=[abs(point.omega),freqs];
        end
        if ~isempty(freqs)
            for i=1:length(freqs)
                [i1,i2]=min(abs(real(l1))+abs(abs(imag(l1)-freqs(i)))); %#ok<ASGLU>
                l1(i2)=-1;
                [i1,i2]=min(abs(real(l1))+abs(abs(imag(l1)+freqs(i)))); %#ok<ASGLU>
                l1(i2)=-1;
            end
        end
        % look for non-real complex pair closest to imaginary axis
        l1=l1(imag(l1)>0);
        if isempty(l1)
            error('P_TOHOPF: no good pair of complex roots found.');
        end
        [i1,i2]=min(abs(real(l1))); %#ok<ASGLU>
        omega=imag(l1(i2));
        x=point.x;
    case 'psol'
        x=sum(point.profile,2)/size(point.profile,2);
        omega=2*pi/point.period;
    otherwise
        try
            conversion=str2func([point.kind,'_tohopf']);
            hopf=conversion(funcs,point,freqs);
            return
        catch ME
            error('p_tohopf: point type %s not supported',point.kind);
        end
end

hopf.x=x;
D=ch_matrix(funcs,x,point.parameter,1i*omega);
[E1,E2]=eig(D);
[i1,i2]=min(abs(diag(E2))); %#ok<ASGLU>
hopf.v=E1(:,i2);
hopf.omega=omega;
end

