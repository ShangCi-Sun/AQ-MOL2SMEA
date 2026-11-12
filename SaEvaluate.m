% Author: Shangci Sun
% Email: sunshangci_stu@163.com
function objs = SaEvaluate(subproblemList,decs)



    [N,nLinear] = size(decs);
    objs = zeros(N,nLinear);
    for k = 1 : nLinear
        objs(:,k) = subproblemList{k}.fobj(decs(:,k));
    end
end