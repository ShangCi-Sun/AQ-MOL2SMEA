% Author: Shangci Sun
% Email: sunshangci_stu@163.com
function [Archive,currFEs,tArc,Population,MaxObj,MinObj] = UpdateArchiveCLT(Archive,Problem,popNew,remain,currFEs,tArc,Population,MaxObj,MinObj)



    newN = size(popNew,1);
    D    = Problem.D;
    BU   = Problem.upper;
    BD   = Problem.lower;
    
    %% Update Archive
    if newN == 0 || remain <= 0
        return;
    elseif newN <= remain
        % Evaluate all candidate solutions with true function
        Dec        = popNew(:,1:D).*repmat(BU-BD,newN,1) + repmat(BD,newN,1);
        NEWOFF     = Problem.Evaluation(Dec);
        Population = [Population,NEWOFF];
        DECS       = (NEWOFF.decs-repmat(Problem.lower,length(NEWOFF),1))./repmat(Problem.upper-Problem.lower,length(NEWOFF),1);
        MinObj     = min([MinObj;NEWOFF.objs],[],1);
        MaxObj     = max([MaxObj;NEWOFF.objs],[],1);
        tArc       = [tArc;DECS, NEWOFF.objs];
        Archive    = [Archive;DECS,NEWOFF.objs];
        currFEs    = currFEs + newN;
    else
        % Evaluate top 'remain' candidate solutions with true function
        [~,idx]    = sort(popNew(:,D+1),'ascend');
        topRemain  = idx(1:remain);
        Dec        = popNew(topRemain,1:D).*repmat(BU-BD,remain,1) + repmat(BD,remain,1);
        NEWOFF     = Problem.Evaluation(Dec);
        Population = [Population,NEWOFF];
        DECS       = (NEWOFF.decs-repmat(Problem.lower,length(NEWOFF),1))./repmat(Problem.upper-Problem.lower,length(NEWOFF),1);
        MinObj     = min([MinObj;NEWOFF.objs],[],1);
        MaxObj     = max([MaxObj;NEWOFF.objs],[],1);
        tArc       = [tArc;DECS, NEWOFF.objs];
        Archive    = [Archive;DECS,NEWOFF.objs];
        currFEs    = currFEs + remain;
    end
end