%% Dissmimilarity Function
function D = mIDENT_dissimilarity(X1,X2,Y1,Y2,S1,S2)

%% Filter bands
% The algorithm compares pairwisely the closeset peaks of the two spectra.
% In this context, a maximum gap size of 20 cm-1 is preset.
allX = pdist2(X1,X2); % Calculate pairwise distance of peak positions.
allX(allX > 20) = NaN; % Delete all pairs with a distance greater that 20.

hlpIDX = sum(~isnan(allX),1) > 1; % Check reference for multipletts;
    allX(allX~=min(allX,[],1)) = nan;

hlpIDX = sum(~isnan(allX),2) > 1; % Check sample for multipletts;
    allX(allX~=min(allX,[],2)) = nan;

hlpIDX1 = nansum(allX,2)==0; % Check sample for unused bands.
hlpIDX2 = nansum(allX,1)==0; % Check reference for unused bands.
    
X1(hlpIDX1) = []; % Delete all unused sample peak positions.
X2(hlpIDX2) = []; % Delete all unused reference peak positions.
Y1(hlpIDX1) = []; % Delete all unused sample peak areas.
Y2(hlpIDX2) = []; % Delete all unused reference peak areas.
S1(hlpIDX1) = []; % Delete all unused sample peak significance values.
      S2all = S2; % Save reference significance vector.
S2(hlpIDX2) = []; % Delete all unused reference peak significance values.

%% Significance Weighting
% The Weighting factor W is defined as the cumulated significance S2 of the
% corresponding peak areas Y2 that can be found in the sample spectrum Y1 as well. W = (sum(S2)/sum(S2all));
W = sum(S2)^2/sum(S2all)^2;


%% Check For Minimum Number Of Peaks
% The dissimilarity function only supports spectra with three or more peaks.
% If this is not the case, the dissimilarity is set to maximum. 
if numel(X1) > 2 && numel(X2) > 2

%% Normalise Peaks
% In this section, all common sample and reference peaks are normalised.
% For that reason, all possible and unique peak pairs are combined using
% the Canberra Distance. 
rY1 = pdist(Y1,@(u,v) (u-v)./(u+v))'; % Normalise sample peaks
rY2 = pdist(Y2,@(u,v) (u-v)./(u+v))'; % Normalise reference peaks
rX1a = pdist(X1,@(u,v) u)'; % Create sample peak position vector 1 (normalised).
rX1b = pdist(X1,@(u,v) v)';% Create sample peak position vector 2 (normalised).
rX2a = pdist(X2,@(u,v) u)'; % Create reference peak position vector 1 (normalised).
rX2b = pdist(X2,@(u,v) v)';% Create reference peak position vector 2 (normalised).
rS1a = pdist(S1,@(u,v) u)'; % Create sample peak significance vector 1 (normalised).
rS1b = pdist(S1,@(u,v) v)';% Create sample peak significance vector 2 (normalised).
rS2a = pdist(S2,@(u,v) u)'; % Create reference peak significance vector 1 (normalised).
rS2b = pdist(S2,@(u,v) v)';% Create reference peak significance vector 2 (normalised).
rS1 = min([rS1a rS1b],[],2); % Calculate sample peak significance vector (normalised).
rS2 = min([rS2a rS2b],[],2); % Calculate reference peak significance vector (normalised).

%% Calculate distances between signals
% The dinstance or dissimilarity includes peak position and area
% differences. In addidtion, a weighting vector is used, which corresponds
% to peak signifiances.
dx = sqrt((rX1a-rX2a).^2 + (rX1b-rX2b).^2); % Normalised peak position difference.
dy = sqrt((rY1-rY2).^2); % Normalised peak area difference.
dw = mean([rS1 rS2],2); % Weighting vector

k1 = 1.2; k2 = 9; k3 = 4.7; k4 = 5.3; % Empirical factors
dz = 1 - 1./((exp(k1*dx-k2)+1).*(exp(k3*dy-k4)+1)); % Normalised peak distance.
D = sum(dz.*dw/sum(dw))^W; % Calculate weighted dissimilarity.
else % If sample and reference have no common vibriational band
    D = 1;
end
end % closing function