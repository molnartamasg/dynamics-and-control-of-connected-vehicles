function result_array=extract_from_tr(trPO_array,component)
%% extract components from extended solution branch or point array
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
%% check if input is branch rather than point array
if ~isfield(trPO_array,'kind') && isfield(trPO_array,'point')
    trPO_array=trPO_array.point;
end
dim=size(trPO_array(1).profile,1)/3;
type={'kind','solution','eigenvector','omega'};
for i=1:length(trPO_array)
    trPO=trPO_array(i);
    switch component
        case type{1} % kind
            result_array='PDorTorusBif';
            break
        case type{2} % solution
            result=trPO;
            result.profile=result.profile(1:dim,:);
            result.parameter=result.parameter(1:end-2);
        case type{3} % eigenvector
            result=trPO;
            result.profile=result.profile(dim+1:end,:);
            result.parameter=result.parameter(end-1:end);
        case type{4} % omega
            result=trPO.parameter(end-1);
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