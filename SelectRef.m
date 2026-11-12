% Author: Shangci Sun
% Email: sunshangci_stu@163.com
function refPoints = SelectRef(nLinear,xmean,sigma,B,Diag)
% This function is used to select several pair of reference points


    
    D = length(xmean);
    % Select reference points
    refPoints = zeros(2*nLinear,D);
    for k = 1 : 2*nLinear
        refPoints(k,:) = (xmean + sigma * B * (Diag.*randn(D,1)))';
    end
    % Repair invalidate reference points, default BU and BD is 1,0
    refPoints(refPoints<0) = 1e-6;
    refPoints(refPoints>1) = 1-(1e-6);
end