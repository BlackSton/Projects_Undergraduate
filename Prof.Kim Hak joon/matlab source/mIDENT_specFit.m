%% Automated Data Preprocessing
function z = mIDENT_specFit(x,y)

%% Baseline Correction
% In a first step, baseline is corrected by differentiation using the 
% Savitzky--Golay--Algorithm. The polynom order "K" is set to K=2 and
% the window size "F" is set to F=2.
K = 2; % polynom order
F = 7; % windowsize
[~, g] = sgolay(K, F); % Savitzky-Golay-Polynom
    g = g(:, 2);
dx = mean(diff(x)); % step size
N = size(x, 1); % number of data points
HalfWin  = ((F+1)/2) - 1;
SG1 = zeros(N-(F+1)/2, size(x,2));

for n = HalfWin+1:N-HalfWin-1
  SG1(n, :) = dot(g, y(n - HalfWin:n + HalfWin, :));
end

yDeriv = SG1./dx; % 1. order derivative of y 
xDeriv = x(1:size(y, 1)-HalfWin-1, :); % corresponding x vector


%% Noise Level Characterisation
% In order to characterise the noise level of the differentiated spectrum,
% a histogram analysis is performed. Considering that, the distribution of
% all differentiated absorbances is fitted using a Lorentzian Model. We
% define the level of noise based on the peak width, which encloses 68,3%
% of the entire peak area.
binwidth = 3.5*std(yDeriv)/size(yDeriv,1)^(1/3); % bin width, based on Scott's rule
nbins = abs(max(yDeriv)-min(yDeriv))/binwidth; % corresponding number of bins
xx = linspace(min(yDeriv),max(yDeriv),nbins); % bin vector (x data)
yy = interp1(unique(yDeriv),1:numel(unique(yDeriv)),xx,'linear'); % bin vector (y data) 

%%
% Lorentzian Model
myfunc =  @(c,a,s,w,x) a./(((x-c).^2.*(exp(-s*(x-c)) + 1).^2)./w.^2 + 1);

%%
% Hitherto, bin vectors describe the cummultative distribution. To apply
% a Lorentzian Model, those data have to be differentiated. Additionally,
% the edges, which contain vibrational band information have to be cutted
% off, considering that this procedure characterises noise level.
xx2 = xx(1:end-1)';
yy2 = diff(yy)'; % differentiate bin vector

hlp = cumsum(yy2)./sum(yy2); % normalised help vector
yy2 = yy2(hlp>.05 & hlp<.95); % cut bin vector (y data)
xx2 = xx2(hlp>.05 & hlp<.95); % cut bin vector (x data)

%%
% Fit the noise distribution with the predefined Lorentzian Model.
myfit = fit(xx2,yy2,myfunc,...
            'StartPoint',[0 max(yy2) 0 0],...
            'Upper',[.1 max(yy2)*1.2 1 1],...
            'Lower',[-.1 max(yy2)*.8 -1 0],...
            'Algorithm','Trust-region');
noise = myfit.w*tan(0.34135*pi); % As mentioned, noise is defined by peak width.


%% Extract Fingerprint
% Assuming that most characteristic vibrational bands can be obtained from
% fingerprint region, this part of the spectrum is extracted. For further data
% evaluation limits of 700cm-1 to 1900cm-1 are predefined.
lmts = [700 1900];
yDeriv = yDeriv(xDeriv>lmts(1) & xDeriv<lmts(2)); % extract fingerprint: y data
xDeriv = xDeriv(xDeriv>lmts(1) & xDeriv<lmts(2)); % extract fingerprint: x data


%% Detect Peaks
% In a first step of peak detection, local extrema have to be determined
% using the 'findpeaks' function. In order to filter some very small peaks,
% a minimum peak height ('noise' level) is defined as threshold.
[iMAX,locMAX] = findpeaks(yDeriv,xDeriv);%,'minpeakheight',0*noise); % find local maxima
[iMIN,locMIN] = findpeaks(-yDeriv,xDeriv);%,'minpeakheight',0*noise); % find local minima

%%
% In a next step, every local maximum is connected to a local minimum in
% its close proximity. These pairs will define the inflection points of the
% vibrational bands in the corresponding raw spectrum.
[x0,idx] = sortrows([locMAX ; locMIN ],1); % sorted vector of all extreme value positions
intensity = [iMAX ones(size(iMAX,1),1) ; iMIN zeros(size(iMIN,1),1)]; % vector of all extreme values
intensity = intensity(idx,:); % sorted vector of all extreme values, based on their positions
opt = [x0(:,1) intensity]; % sorted matrix, based on extreme value position
                           % [position, intensity, type]
                           % type:1 --> local maximum
                           % type:0 --> local minimum
opt_clean = [];
%%
% The peak matrix 'opt' contains some 'noisy' local extrema that do not
% describe an inflection point of a significant vibrational band. Those
% signals have to be eliminated in the following loop.
while ~isempty(opt)
    if opt(1,3) == 1
        hlp = find(opt(:,3)-1);
    else
        hlp = find(opt(:,3));
    end
    if ~isempty(hlp)
        [~,idxOpt] = max(opt(1:hlp-1,2));
        opt_clean = [opt_clean; opt(idxOpt,:)];
        opt(1:hlp-1,:) = [];
    else
        [~,idxOpt] = max(opt(1:end,2));
        opt_clean = [opt_clean; opt(idxOpt,:)];
        opt = [];
    end
end

%%
% In absorbance mode, every vibrational band starts with a local maximum
% followed by a local minimum. Considering that, it has to be checked if
% the peak matrix 'opt' starts with a local minimum or ends with a local
% maximum, respectively.
if opt_clean(1,3) == 0
    opt_clean(1,:) = [];
end
if opt_clean(end,3) == 1
    opt_clean(end,:) = [];
end

%%
% Inflection points are combined to new matrices: [position of maximum,
% position of minimum] & [intensity of maximum, intensity of minimum]
inflPtsPos = [opt_clean(opt_clean(:,3)==1,1) opt_clean(opt_clean(:,3)==0,1)];
inflPtsInt = [opt_clean(opt_clean(:,3)==1,2) opt_clean(opt_clean(:,3)==0,2)];

%% Initial Fit Parameters
% In order to perform curve fitting using first order derivative of
% asymmetric pseudo Voigt function, suitable initial parameters have to be
% estimated.

%%
% Vibtrational band position: x0
x0 = zeros(size(inflPtsPos,1),1);
for u = 1:size(inflPtsPos,1)
    tmpX = xDeriv(xDeriv>=inflPtsPos(u,1) & xDeriv<=inflPtsPos(u,2));
    [~,idx] = min(abs(yDeriv(xDeriv>=inflPtsPos(u,1) & xDeriv<=inflPtsPos(u,2))));
    x0(u) = tmpX(idx);
end
%%
% Vibrational band width: w
w = (inflPtsPos(:,2)-inflPtsPos(:,1))/2;

%%
% Vibrational band height: a
a = max(inflPtsInt*2.*w,[],2);

%%
% Asymmetry factor: s
R = (inflPtsInt(:,1)-inflPtsInt(:,2))./(inflPtsInt(:,1)+inflPtsInt(:,2)); % Camberra Distance of inflection point intensities 
    a1 = 2.007663533725664; % empirical fit parameter
    b1 = -56.245235357193230; % empirical fit parameter
    c1 = 0.030809453318739; % empirical fit parameter
s = -log(a1./(R.*exp(c1.*w) + 1) - 1)./b1; % empiricial estimation of f
s = real(s)-imag(s);

%%
% Shape factor: f
f = x0*0+.2; % all shape factors are preset to 0.2

%%
% Combine all initial parameters in one matrix
initial = [x0 a w f s];

%% 
% Filter peaks based on minimum and maximum peak width and minimum peak height
% minheight = min([noise*4,0.0023]); original function
minwidth = 2*mean(diff(xDeriv));
maxwidth = 100;
minheight = noise*3;
initial(...
    initial(:,3)<minwidth |...
    initial(:,3)>maxwidth |...
    initial(:,2)<minheight,:) = [];%|...
    %initial(:,2)./initial(:,3)<.0013,:) = [];

%%
% Set upper and lower limits of fit parameters.
up = [initial(:,1) + minwidth,...
      initial(:,2) * 10,...
      initial(:,3) * 10,...
      ones(size(initial,1),1),...
      ones(size(initial,1),1) * 0.2];
lo = [initial(:,1) - minwidth,...
      initial(:,2) * .1,...
      initial(:,3) * .1,...
      zeros(size(initial,1),1),...
      ones(size(initial,1),1) * -0.2];
%%
% Reshape parameter matrices to parameter row vectors.
initial = reshape(initial',1,[]);
up = reshape(up',1,[]);
lo = reshape(lo',1,[]);

%% Curve Fitting
% In a next step, the incomming raw data (x,y) are fitted using a non-linear
% regression. In this context, first order derivative of asymmetric pseudo
% Voigt function is chosen as fit function.

%%
% Considering that raw data (x,y) contains not a single but multiple vibrational
% bands, a cumulative fit function has to be defined.
n = length(initial)/5; % number of peaks
FCNprefix = 'ffunc9(x,';
idxArray = arrayfun(@(u) ['p' sprintf('%03d',u)],1:n,'uni',0);
parArray = {'_p1_c';'_p2_A';'_p3_w';'_p4_mu';'_p5_a'};
varArray = cell(5,n);
for u = 1:n
    for v = 1:5
        varArray{v,u} = [idxArray{u} parArray{v} ','];
    end
end
FCNstr = strjoin(reshape(varArray,1,[]));
FCNstr(end) = [];
FCNstr = [FCNprefix FCNstr ')']; % fit function
%%
% Options for curve fitting:
ft = fittype(FCNstr);
options = fitoptions(ft);
    options.StartPoint = initial;
    options.Lower = lo;
    options.Upper = up;
    options.Robust = 'Bisquare';
    options.Display = 'off';
%%
% Fit
[FIT.curve, FIT.gof] = fit(xDeriv,yDeriv,ft,options);
%% Output: Positions, Areas and Weightings
FITcoefficients = coeffvalues(FIT.curve); % read fit parameters
heights = FITcoefficients(2:5:end); % fitted peak heighths
widths = FITcoefficients(3:5:end); % fitted peak widths
shapes = FITcoefficients(4:5:end); % fitted peak shapes
asymm = FITcoefficients(5:5:end); % fitted asymmetry factors
%%
% Fitted peak positions 
z.Positions = FITcoefficients(1:5:end)' 

%%
% fitted peak areas, calculated by numerical integration
z.Areas = arrayfun(@(u) ... 
    trapz(...
        linspace(... % x data for integration
            z.Positions(u)-8*widths(u),...
            z.Positions(u)+8*widths(u),...
            500),...
        mypeakfunction2(... % y data for integration
            linspace(...
                z.Positions(u)-8*widths(u),...
                z.Positions(u)+8*widths(u),...
                500),...
            z.Positions(u),...
            heights(u),...
            widths(u),...
            shapes(u),...
            asymm(u)...
            )...
        ),...
    1:n)';

%%
% Fitted peak weightings
%   z.Weightings(i) = z.Areas(i)/sqrt(sum((FIT.data(tmpRange)-yDeriv(tmpRange)).^2)/sum((yDeriv(tmpRange)-mean(yDeriv(tmpRange))).^2)); 
z.Weightings = zeros(n,1);
FIT.data = ffunc9(xDeriv,FITcoefficients); % fitted y data
for i = 1:n
    tmpRange = xDeriv > z.Positions(i) - 1.64*widths(i) &...
               xDeriv < z.Positions(i) + 1.64*widths(i); % the factor 1.64 encloses 90% of peak area
    z.Weightings(i) = z.Areas(i)/sqrt(sum((FIT.data(tmpRange)-y1(tmpRange)).^2)/sum((y1(tmpRange)-mean(y1(tmpRange))).^2)); 
end

end % closing function