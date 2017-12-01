function runMiia
clc
close all

% choose data
%----------------------------
% dataSource = 'tutorial';  
dataSource = 'friedman';  
% dataSource = 'glv';
%----------------------------

switch dataSource
    case {'tutorial','friedman'}
        fileName=[pwd,'\data\smallData\tutorialFriedman.xlsx'];
        miia=goMiia(dataSource,fileName);
        aB_pred=miia.aB; aC_pred=miia.aC;
        aB_pred
        aC_pred
        disp('No true values of interaction coefficients for comparison are available.')
    case 'glv'
        %-- don't change
        r=1; % growth rate
        alpha=0.1; % std of binary interaction coefficients 
        %-- choose specific values for S, beta, and trial
        S=15; % species number: 3:15
        beta=0.5; % std of complex interaction coefficients: [0:0.02:0.1,0.15,0.2,0.25,0.3:0.1:1]
        trial=1; % # of trial: 1:30
        fileName=[pwd,'\data\glvData\','glv_',num2str(S),'_',num2str(r),'_',num2str(alpha),'_',num2str(beta),'_',num2str(trial),'.mat'];
        
        miia=goMiia(dataSource,fileName);

        % true interaction coefficients set in the glv model
        load(fileName);
        aB_true=model.aB;
        aC_true=model.aC;

        % interaction coefficients predicted by MIIA
        aB_pred=miia.aB; aC_pred=miia.aC;

        % print results
        aB_true
        aC_true
        disp('---')
        aB_pred
        aC_pred

        % correlation, similarity, and figures 
        postanalysis(miia,model)
end

%--------------------------------------------------------------------------
function miia=goMiia(dataSource,fileName)

while 1
    [S,sName,xA,xB_all,xC_all,A,Ap] = readBasicData(dataSource,fileName);
    % check solution range
    xMax = max([max(xA(:)),max(xB_all(:)),max(xC_all(:))]);
    xMin = min([min(xA(:)),min(xB_all(:)),min(xC_all(:))]);
    if xMax<=1.e5 && xMin>=-eps        
        break
    else
        disp([xMin,xMax])
        disp('check the range of data')
        pause
    end
end

aB=getBinaryCoeff(S,xA,xB_all);
miia.aB=aB;

nC = size(xC_all,1);
for i=1:nC
    aC=NaN(S,S);
    xC=xC_all(i,:);
    iiComplex = find(~isnan(xC));
    xA_ = diag(xA(iiComplex,iiComplex));
    xC_ = xC(iiComplex);
    aB_ = aB(iiComplex,iiComplex);
    aC_ = getComplexCoeff(aB_,xA_,xC_);
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

% default species name
if isempty(sName)
    for i=1:S
        sName_{i}=['x',num2str(i)]; 
    end
    sName=cell(sName_);
end

%--------------------------------------------------------------------------
function aB = getBinaryCoeff(S,xA,xB_all)
nB = nchoosek(S,2);
aB = zeros(S,S);

for j=1:nB
    iiB = find(~isnan(xB_all(j,:)));
    delx(1,1) = xB_all(j,iiB(1))-xA(iiB(1),iiB(1));
    if xB_all(j,iiB(2))>1e-6
        aB(iiB(1),iiB(2)) = delx(1)/xB_all(j,iiB(2));
    end
    delx(2,1) = xB_all(j,iiB(2))-xA(iiB(2),iiB(2));
    if xB_all(j,iiB(1))>1e-6
        aB(iiB(2),iiB(1)) = delx(2)/xB_all(j,iiB(1));
    end
end

%--------------------------------------------------------------------------
function aC = getComplexCoeff(aB,xA,xC)
S = length(xC);
aC = zeros(S,S);

for i=1:S
    delx = -(xC(i)-xA(i));
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

%--------------------------------------------------------------------------
function postanalysis(miia,model)
S=length(model.R);

% true coefficients
aB=model.aB; aC=model.aC;

% predicted coefficients
aBhat=miia.aB; aChat=miia.aC;

% weighted by populations
xC=model.xAll(end,:);
aChatxC=aChat.*repmat(xC,S,1);
aCxC=aC.*repmat(xC,S,1);

% make column vectors
aB_=aB(:);aB_(1:S+1:end)=[];
aBhat_=aBhat(:);aBhat_(1:S+1:end)=[];
aChat_=aChat(:);aChat_(1:S+1:end)=[];
aC_=aC(:);aC_(1:S+1:end)=[];
aChatxC_=aChatxC(:);aChatxC_(1:S+1:end)=[];
aCxC_=aCxC(:);aCxC_(1:S+1:end)=[];

% cosine similarity: unweigted vs. population weighted  
cosineBC=(aB_'*aC_)/norm(aB_)/norm(aC_);
cosineBBhat=(aB_'*aBhat_)/norm(aB_)/norm(aBhat_);
cosineCChat=(aC_'*aChat_)/norm(aC_)/norm(aChat_);
cosineCChatxC=(aCxC_'*aChatxC_)/norm(aCxC_)/norm(aChatxC_);
cosine=[cosineBC,cosineBBhat,cosineCChat,cosineCChatxC];

% plot the results
% -- normalize
aB_=aB_/max(abs(aB_));
aBhat_=aBhat_/max(abs(aBhat_));
aC_=aC_/max(abs(aC_));
aChat_=aChat_/max(abs(aChat_));
aCxC_=aCxC_/max(abs(aCxC_));
aChatxC_=aChatxC_/max(abs(aChatxC_));
%-- plot
figure
subplot(2,2,1)
plot(aB_,aC_,'o')
xlabel('a_{B,true}'); ylabel('a_{C,true}')
title(['Cos similarity = ',num2str(cosine(1))])
axis square
subplot(2,2,2)
plot(aB_,aBhat_,'o')
xlabel('a_{B,true}'); ylabel('a_{B,miia}')
title(['Cos similarity = ',num2str(cosine(2))])
axis square
subplot(2,2,3)
plot(aC_,aChat_,'o')
xlabel('a_{C,true}'); ylabel('a_{C,miia}')
title(['Cos similarity = ',num2str(cosine(3))])
axis square
subplot(2,2,4)
plot(aCxC_,aChatxC_,'o')
xlabel('a_{C,true} (weighted)'); ylabel('a_{C,miia} (weighted)')
title(['Cos similarity = ',num2str(cosine(4))])
axis square
