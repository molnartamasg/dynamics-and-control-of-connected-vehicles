function result_array=extract_from_POfold(pfold_array,component,npar)
%% extract components from pfold solution branch or point array
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
%% check if input is branch rather than point array
if ~isfield(pfold_array,'kind') && isfield(pfold_array,'point')
    pfold_array=pfold_array.point;
end
%% extract named components
dim=size(pfold_array(1).profile,1)/2;
type={'kind','solution','nullvector','delays'};
for i=1:length(pfold_array)
    pfold=pfold_array(i);
    switch component
        case type{1} %'kind'
            result_array='POfold';
            break
        case type{2} %'solution'
            result=pfold;
            result.profile=result.profile(1:dim,:);
            result.parameter=result.parameter(1:npar);
        case type{3} %'nullvector'
            result=pfold;
            result.profile=result.profile(dim+1:end,:);
            result.parameter=result.parameter(npar+1);
        case type{4} %'delays'
            result.values=pfold.parameter(npar+3:end);
        otherwise
            fprintf('known component types:\n');
            for k=1:length(type)
                fprintf('%s\n',type{k});
            end
            result_array=[];
            break;
    end
    result_array(i)=result; %#ok<AGROW>
end
end
