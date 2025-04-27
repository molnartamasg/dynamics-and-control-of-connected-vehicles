function [options,passed_on]=dde_set_options(defaults,userargs,pass_on)
%% parses varargin and assigns fields of structure options
% unknown arguments are passed on into cell array passed_on if pass_on is
% present and non-empty or false, otherwise and error message is generated
%
% (c) DDE-BIFTOOL v. 3.1.1(19), 11/04/2014
%
passed_on={};
% wrap cell arguments to avoid generating multiple structs
if isstruct(defaults)
    options=defaults;
elseif iscell(defaults)
    for i=1:length(defaults)
        if iscell(defaults{i})
            defaults{i}=defaults(i);
        end
    end
    options=struct(defaults{:});
else
    error('defaults not recognized\n');
end
if nargin<3 || isempty(pass_on)
    pass_on=false;
end
for i=1:2:length(userargs)
    if isfield(options,userargs{i})
        options.(userargs{i})=userargs{i+1};
    else
        if ~pass_on
            error('option ''%s'' not recognized\n',userargs{i});
        else
            passed_on={passed_on{:},userargs{i},userargs{i+1}}; %#ok<CCAT>
        end
    end
end
end
