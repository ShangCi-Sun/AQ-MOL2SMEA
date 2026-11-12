% Author: Shangci Sun
% Email: sunshangci_stu@163.com
function candidates = AQ_MOEAOptimizeOS(subproblemList,popSize,maxIter)
% Optimize each subproblem by multi-objective algorithm



    %% Get dimension and boundary
    D = numel(subproblemList);
    [BU,BD] = deal(zeros(1,D));
    for k = 1 : D
        BU(k) = subproblemList{k}.ub;
        BD(k) = subproblemList{k}.lb;
    end
    
    %% Generate the random population
    population.decs = rand(popSize,D).*repmat(BU-BD,popSize,1) + repmat(BD,popSize,1);
    population.objs = SaEvaluateOS(subproblemList,population.decs);
    
    %% Optimize by multi-objective algorithm
    candidates = NSGAIIIopt(subproblemList,BU,BD,population,popSize,maxIter);
    [~,stds]   = SaEvaluateOS(subproblemList,candidates(:,1:D));
    candidates = [candidates, stds];
end