% Author: Shangci Sun
% Email: sunshangci_stu@163.com
function [objs, stds] = SaEvaluateOS(subproblemList,decs)



    [N,nLinear] = size(decs);
    objs = zeros(N,nLinear);
    stds = zeros(N,nLinear);
    for k = 1 : nLinear
        [objs(:,k), stds(:,k)] = subproblemList{k}.fobj(decs(:,k));
    end
end