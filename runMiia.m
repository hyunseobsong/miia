function runMiia2
clc
% close all
%==========================================================================
% Original MIIA: when species abundance data in axenic, binary, and complex
communities is fully available.
% inferenceRule = 'miia1';
% Scaled MIIA: when only data from axenic cultures are unavailable, to pre-
dict neighbor-dependent interactions in a relative sense
inferenceRule = 'miia2';
%------------------------
dataSource = 'kato';
% dataSource = 'wang';
% dataSource = 'friedman';
% dataSource = 'tutorial';
% dataSource = 'glv';
%------------------------
force2Zero = true; % if set unidentifiable binary coeffs to 0
%==========================================================================

if dataSource=="glv"
    %-- don't change
    r=1; % growth rate
    alpha=0.1; % std of binary interaction coefficients 
    %-- choose specific values for S, beta, and trial
    S=15; % species number: 3:15
    beta=0.5; % std of complex interaction coefficients: [0:0.02:0.1,0.15,0.2,0.25,0.3:0.1:1]
    trial=1; % # of trial: 1:30
    fileName=['data/glvData/','glv_',num2str(S),'_',num2str(r),'_',num2str(alpha),'_',num2str(beta),'_',num2str(trial),'.mat'];
else
    fileName='data/demodata.xlsx';
end

miia=goMiia(dataSource,inferenceRule,force2Zero,fileName);
miia.aB
miia.aC

save([dataSource,'.',inferenceRule,'.mat'],'miia')

%--------------------------------------------------------------------------
function miia=goMiia(dataSource,inferenceRule,force2Zero,fileName)

while 1
    [S,sName,xA,xB_all,xC_all,A,Ap] = readBasicData(dataSource,fileName);
    % check solution range
    xMax = max([max(xA(:)),max(xB_all(:)),max(xC_all(:))]);
    xMin = min([min(xA(:)),min(xB_all(:)),min(xC_all(:))]);
    if xMax<=1.e10 && xMin>=-eps        
        break
    else
        disp([xMin,xMax])
        disp('check the range of data')
        pause
    end
end

aB=getBinaryCoeff(inferenceRule,force2Zero,S,xA,xB_all);
miia.aB=aB;

nC = size(xC_all,1);
for i=1:nC
    aC=NaN(S,S);
    xC=xC_all(i,:);
    iiComplex = find(~isnan(xC));
    xA_ = diag(xA(iiComplex,iiComplex));
    xC_ = xC(iiComplex);
    aB_ = aB(iiComplex,iiComplex);
    aC_ = getComplexCoeff(inferenceRule,aB_,xA_,xC_);
    aC(iiComplex,iiComplex)=aC_;
    miia.aC(1:S,1:S,i)=aC;
end


%--------------------------------------------------------------------------
function [S,sName,xA,xB_all,xC_all,A,Ap] = readBasicData(dataSource,fileName)

A = [];
Ap = [];

sName=[];
switch dataSource
    case 'glv'
        load(fileName); data=model.xAll;
    otherwise
        [data,sName] = xlsread(fileName,dataSource);
end

[nExp,S] = size(data); % number of species
nB = nchoosek(S,2);
% get x_i in axenic cultures 
for i=1:nchoosek(S,1)
    xA(i,:)=data(i,:);
end
% get x_i in binary cultures
for i=1:nB
    xB_all(i,:)=data(S+i,:);
end
% get x_i in complex cultures
for i=1:nExp-(S+nB)
    xC_all(i,:)=data(S+nB+i,:);
end

% default name
if isempty(sName)
    for i=1:S
        sName_{i}=['x',num2str(i)]; 
    end
    sName=cell(sName_);
end

%--------------------------------------------------------------------------
function aB = getBinaryCoeff(inferenceRule,force2Zero,S,xA,xB_all)
nB = nchoosek(S,2);
aB = NaN(S,S);

for j=1:nB
    iiB = find(~isnan(xB_all(j,:)));
    delx(1,1) = xB_all(j,iiB(1))-xA(iiB(1),iiB(1));
    
    if force2Zero
        aB(iiB(1),iiB(2)) = 0;
        aB(iiB(2),iiB(1)) = 0;
    end
    
    if xB_all(j,iiB(2))>1e-6
        if inferenceRule=="miia1"
            if xA(iiB(1),iiB(1))>1e-6
                aB(iiB(1),iiB(2)) = delx(1)/xB_all(j,iiB(2))/xA(iiB(1),iiB(1));
            end
        elseif inferenceRule=="miia2"
            aB(iiB(1),iiB(2)) = delx(1)/xB_all(j,iiB(2));
        end
    end
    delx(2,1) = xB_all(j,iiB(2))-xA(iiB(2),iiB(2));
    if xB_all(j,iiB(1))>1e-6
        if inferenceRule=="miia1"
            if xA(iiB(2),iiB(2))>1e-6
                aB(iiB(2),iiB(1)) = delx(2)/xB_all(j,iiB(1))/xA(iiB(2),iiB(2));
            end
        elseif inferenceRule=="miia2"
            aB(iiB(2),iiB(1)) = delx(2)/xB_all(j,iiB(1));
        end
    end
end

%--------------------------------------------------------------------------
function aC = getComplexCoeff(inferenceRule,aB,xA,xC)
S = length(xC);
aC = NaN(S,S);

for i=1:S
    if inferenceRule=="miia1"
        delx = -(xC(i)-xA(i))/xA(i);
    elseif inferenceRule=="miia2"
        delx = -(xC(i)-xA(i));
    end
    xC_ = xC; xC_(i)=[];
    aB_ = aB(i,:); aB_(i)=[];
    xC_ = xC_(:); aB_ = aB_(:); % column vectors
    num_ = xC_'*aB_ + delx;
    den_ = sumsqr(xC_);
    common = num_/den_;
    for j=1:S-1
        aC_(1,j) = aB_(j)-xC_(j)*common;
    end
    iiC_ = 1:S; iiC_(i) = [];
    aC(i,iiC_)=aC_;
end
