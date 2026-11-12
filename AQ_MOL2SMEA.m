% Author: Shangci Sun
% Email: sunshangci_stu@163.com
classdef AQ_MOL2SMEA < ALGORITHM
% <2025> <multi> <real> <expensive> <large/none>
% Multi-objective linear subspace surrogate modeling assisted evolutionary algorithm
% NLinear --- 8 --- Number of linear subspaces



    methods
        function main(Algorithm,Problem)
            %% Generate random population
            NLinear    = Algorithm.ParameterSet(8);
            NT         = 2*Problem.D;
            MaxFEs     = Problem.maxFE - NT;
            Dec        = lhsdesign(NT,Problem.D).*repmat(Problem.upper-Problem.lower,NT,1) + repmat(Problem.lower,NT,1);
            Population = Problem.Evaluation(Dec);
            Objs       = Population.objs;
            MaxObj     = max(Objs,[],1);
            MinObj     = min(Objs,[],1);
            TArchive   = [Dec,Population.objs];
            TArchive(:,1:Problem.D) = (TArchive(:,1:Problem.D)-repmat(Problem.lower,NT,1))./repmat(Problem.upper-Problem.lower,NT,1);
            
            %% Set parameters
            [~,MaxIter] = Algorithm.ParameterSet(40,30);
            BasePop = Problem.N;
            PopSize = BasePop;
            % Parameter of surrogate model
            CurrFEs = 0;
            
            %% Cluster number is 5
            CLTn     = 5;
            [W,CLTn] = UniformPoint(CLTn, Problem.M);
            %% Linear subspace quality information

            SubspaceReward  = zeros(CLTn,NLinear);

            SubspaceCount   = zeros(CLTn,NLinear);

            SubspaceSuccess = zeros(CLTn,NLinear);

            SubspaceQuality = zeros(CLTn,NLinear);
            % Parameter in CMA-ES
            XMeanCell    = cell(1, CLTn);
            SigmaCell    = cell(1, CLTn);
            PathConvCell = cell(1, CLTn);
            PathSigCell  = cell(1,CLTn);
            BCell        = cell(1,CLTn);
            DiagCell     = cell(1,CLTn);
            ConvCell     = cell(1,CLTn);
            InvSqrtCCell = cell(1,CLTn);
            NObjs        = (TArchive(:, end-Problem.M+1:end) - repmat(MinObj,size(TArchive,1),1))./repmat(MaxObj-MinObj,size(TArchive,1),1);
            [~,grp]      = min(pdist2(W, NObjs,'cosine'), [], 1);
            for k = 1 : CLTn
               
                if any(grp==k)
                    XMeanCell{k} = mean(TArchive(grp==k,1:Problem.D),1)';
                else
                    XMeanCell{k} = mean(TArchive(:,1:Problem.D),1)';
                end
                SigmaCell{k}    = 0.5;
                PathConvCell{k} = zeros(Problem.D,1);
                PathSigCell{k}  = zeros(Problem.D,1);
                BCell{k}        = eye(Problem.D,Problem.D);
                DiagCell{k}     = ones(Problem.D,1);
                ConvCell{k}     = BCell{k}*diag(DiagCell{k}.^2)*BCell{k}';
                InvSqrtCCell{k} = BCell{k}*diag(DiagCell{k}.^-1)*BCell{k}';
            end
            EigenEval = 0;
            ChiD      = Problem.D^0.5*(1-1/(4*Problem.D)+1/21*Problem.D^2);

            % Selection
            Lambda  = 4*NLinear;
            Mu      = Lambda/2;
            Weights = log(Mu+1/2) - log(1:Mu)';
            Mu      = floor(Mu);
            Weights = Weights/sum(Weights);
            MuEff   = sum(Weights)^2/sum(Weights.^2);

            % Adaptation
            CumConv = (4+MuEff/Problem.D)/(Problem.D+4+2*MuEff/Problem.D);
            CumSig  = (MuEff+2)/(Problem.D+MuEff+5);
            CRate   = 2/((Problem.D+1.3)^2+MuEff);
            CMu     = min(1-CRate,2*(MuEff-2+1/MuEff)/((Problem.D+2)^2+MuEff));
            Damps   = 1 + 2*max(0,sqrt((MuEff-1)/(Problem.D+1))-1);

            %% Optimization
            NextPrintFE = Problem.FE + 500;
            while Algorithm.NotTerminated(Population)
                LoopStart = tic;
                FEStart = Problem.FE;
                TotalConstruct = 0;
                TotalOptimize = 0;
                TotalUpdate = 0;
                TotalSelect = 0;
                TotalCandidate = 0;
                TotalCMA = 0;
                TotalRewardSort = 0;
                Reward = SubspaceReward./(SubspaceCount+1e-6);

                Success = SubspaceSuccess./(SubspaceCount+1e-6);

                Explore = 1./(SubspaceCount+1);


                SubspaceQuality = ...
                    0.5*Reward + ...
                    0.3*Success + ...
                    0.2*Explore;
      

                NObjs   = (TArchive(:, end-Problem.M+1:end) - repmat(MinObj,size(TArchive,1),1))./repmat(MaxObj-MinObj,size(TArchive,1),1);
                [~,grp] = min(pdist2(W, NObjs,'cosine'), [], 1);


                %% Quality-driven adaptive population allocation

                ClusterQuality = mean(SubspaceQuality,2);


                beta = 5;

                ExpQ = exp(beta*(ClusterQuality-max(ClusterQuality)));

                QualityWeight = ExpQ./sum(ExpQ);

                WeightNorm = ...
                    (QualityWeight-min(QualityWeight)) / ...
                    (max(QualityWeight)-min(QualityWeight)+1e-6);
                Pop_min = round(0.5*BasePop);

                Pop_max = round(1.5*BasePop);
                AdaptivePop = round(...
                    Pop_min + ...
                    (Pop_max-Pop_min)*WeightNorm);
               %%
                BaseBudget = 10;

                Budget_min = round(0.5*BaseBudget);

                Budget_max = round(2*BaseBudget);

                AdaptiveFEBudget = round(...
                    Budget_min + ...
                    (Budget_max-Budget_min)*WeightNorm);
               
                %%
                MinPopSafe = 20;

                AdaptivePop = max(AdaptivePop,MinPopSafe);
                AdaptivePop = AdaptivePop + mod(AdaptivePop,2);

                for k = 1 : CLTn
                   if CurrFEs >= MaxFEs || Problem.FE >= Problem.maxFE
                    break;
                end
                    %% Construct several linear space
                    % Select reference points
                    tSelect = tic;
                    RefPoints = SelectRef(NLinear,XMeanCell{k},SigmaCell{k},BCell{k},DiagCell{k});
                    tArc = TArchive(grp==k,:);
                    if sum(grp==k) < 30
                        disth = pdist2(W(k,:),NObjs,'cosine');
                        [~,indx] = sort(disth);
                        tArc = TArchive(indx(1:2*Problem.D),:);
                    end
                    TotalSelect = TotalSelect + toc(tSelect);
                    % Construct process
                    tConstruct = tic;

                    SubproblemList = AQ_ConstructLinearKRGOS( ...
                        RefPoints,tArc,NLinear,MaxObj,MinObj,W(k,:));

                    TotalConstruct = TotalConstruct + toc(tConstruct);

                    %% Optimize process
                    % Optimize by quality-driven adaptive population
                    
                    CurrentPopSize = AdaptivePop(k);

                    tOptimize = tic;

                    Candidates = AQ_MOEAOptimizeOS( ...
                        SubproblemList,CurrentPopSize,MaxIter);

                    TotalOptimize = TotalOptimize + toc(tOptimize);
                    % Select candidate solutions
                    Remain    = MaxFEs - CurrFEs;
                    tCandidate = tic;

                    [Offspring,SourceID] = AQ_SelectCandidateKRGOS( ...
                        SubproblemList,NLinear,Problem.D,Candidates);

                    TotalCandidate = TotalCandidate + toc(tCandidate);

                    Budget = min(AdaptiveFEBudget(k),Remain);

                    if size(Offspring,1) > Budget

                        Offspring = Offspring(1:Budget,:);

                        SourceID  = SourceID(1:Budget);

                    end
                    
                    %% Update Archive

                    tRewardSort = tic;

                    OldFrontNo = NDSort(TArchive(:,Problem.D+1:end),1);
                    OldFrontNum = sum(OldFrontNo==1);

                    TotalRewardSort = TotalRewardSort + toc(tRewardSort);
                    OldArchiveSize = size(TArchive,1); 
                    tUpdate = tic;


                    [TArchive,CurrFEs,tArc,Population,MaxObj,MinObj] = ...
                        UpdateArchiveCLT(TArchive,Problem,Offspring,Remain,...
                        CurrFEs,tArc,Population,MaxObj,MinObj);


                    TotalUpdate = TotalUpdate + toc(tUpdate);

                    tRewardSort = tic;

                    NewObj = TArchive(:,Problem.D+1:end);
                    NewFrontNo = NDSort(NewObj,1);
                    NewFrontNum = sum(NewFrontNo==1);

                    TotalRewardSort = TotalRewardSort + toc(tRewardSort);
                  
                    DeltaReward = max(NewFrontNum-OldFrontNum,0);
                    %% Update subspace quality information

                    for i = 1:length(SourceID)

                        sid = SourceID(i);

                        % number of evaluations
                        SubspaceCount(k,sid)=...
                            SubspaceCount(k,sid)+1;


                        % contribution reward
                        SubspaceReward(k,sid)=...
                            SubspaceReward(k,sid)+DeltaReward/length(SourceID);

                    end
                    %% Calculate offspring success

                    OffspringObj = TArchive(OldArchiveSize+1:end,Problem.D+1:end);

                    assert(size(OffspringObj,1) == numel(SourceID), ...
                        '实际评价解数量与SourceID数量不一致');

                    [OffspringFront,~] = NDSort(OffspringObj,1);


                    for i = 1:length(SourceID)

                        sid = SourceID(i);

                        if OffspringFront(i)==1

                            SubspaceSuccess(k,sid)=...
                                SubspaceSuccess(k,sid)+1;

                        end

                    end


                    % Keep archive ranking for CMA-ES update
                    tCMA = tic;
                    [FrontNo,~] = NDSort(tArc(:,Problem.D+1:end),1);
                    
                    if size(tArc,1) >= Mu
                        in = find(FrontNo==1);
                        id = find(FrontNo~=1);
                        if length(in) >= ceil(Mu/2)
                            if length(id) >= ceil(Mu/2)
                                sni     = randperm(length(in),ceil(Mu/2));
                                sdi     = randperm(length(id),ceil(Mu/2));
                                MeanDec = [tArc(in(sni),1:Problem.D);tArc(id(sdi),1:Problem.D)];
                            else
                                sni     = randperm(length(in),Mu-length(id));
                                MeanDec = [tArc(id,1:Problem.D);tArc(in(sni),1:Problem.D)];
                            end                    
                        else
                            sdi     = randperm(length(id),Mu-length(in));
                            MeanDec = [tArc(in,1:Problem.D);tArc(id(sdi),1:Problem.D)];
                        end
                    else
                        dist    = pdist2(W(k,:), NObjs,'cosine');
                        [~,cdt] = sort(dist,'ascend');
                        MeanDec = TArchive(cdt(1:Mu),1:Problem.D);
                    end

                    XOld         = XMeanCell{k};
                    XMeanCell{k} = MeanDec'*Weights;
                    
                    PathSigCell{k}  = (1-CumSig)*PathSigCell{k} + sqrt(CumSig*(2-CumSig))*MuEff*InvSqrtCCell{k}*(XMeanCell{k}-XOld)/SigmaCell{k};
                    HSig            = norm(PathSigCell{k})/sqrt(1-(1-CumSig)^(2*CurrFEs/Lambda))/ChiD < 1.4+2/(Problem.D+1);
                    PathConvCell{k} = (1-CumConv)*PathConvCell{k} + HSig*sqrt(CumConv*(2-CumConv)*MuEff)*(XMeanCell{k}-XOld)/SigmaCell{k};
                    ArcTmp          = (1/SigmaCell{k})*(MeanDec'-repmat(XOld,1,Mu));
                    
                    ConvCell{k} = (1-CRate-CMu)*ConvCell{k} + CRate*(PathConvCell{k}*PathConvCell{k}'+(1-HSig)*CumConv*(2-CumConv)*ConvCell{k}) + CMu*ArcTmp*diag(Weights)*ArcTmp';
                    
                    SigmaCell{k} = SigmaCell{k}*exp((CumSig/Damps)*(norm(PathSigCell{k})/ChiD-1));
                    if CurrFEs - EigenEval > Lambda/(CRate+CMu)/Problem.D/10
                        EigenEval       = CurrFEs;
                        ConvCell{k}     = triu(ConvCell{k}) + triu(ConvCell{k},1)';
                        [BCell{k},DiagCell{k}] = eig(ConvCell{k});
                        DiagCell{k}     = sqrt(diag(DiagCell{k}));
                        InvSqrtCCell{k} = BCell{k}*diag(DiagCell{k}.^-1)*BCell{k}';
                        if ~isreal(DiagCell{k})
                            DiagCell{k} = real(DiagCell{k});
                        end
                    end
                    TotalCMA = TotalCMA + toc(tCMA);
                end 
                 LoopTime = toc(LoopStart);
                AddedFE = Problem.FE - FEStart;

                % if Problem.FE >= NextPrintFE
                %     fprintf(['FE=%d | Added=%d | Archive=%d | Loop=%.2fs | ',...
                %              'RewardSort=%.2fs | Build=%.2fs | Optimize=%.2fs | ',...
                %              'Update=%.2fs | Select=%.2fs | Candidate=%.2fs | CMA=%.2fs\n'],...
                %         Problem.FE,AddedFE,size(TArchive,1),LoopTime,...
                %         TotalRewardSort,TotalConstruct,TotalOptimize,...
                %         TotalUpdate,TotalSelect,TotalCandidate,TotalCMA);
                % 
                %     NextPrintFE = Problem.FE + 500;
                % end
               
            end
        end
    end
end
 