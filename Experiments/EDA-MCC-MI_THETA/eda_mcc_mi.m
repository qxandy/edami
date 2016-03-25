% =========================================================================
% -------------------------------------------------------------------------
%
%        IMPLEMENTATION OF EDA-MCC BY DONG ET AL.
% -------------------------------------------------------------------------
% =========================================================================
%    Version : 1.0 (7 March 2013)
%    Since   : 27 Feb 2014
%    Authors : Jakramate Bootkrajang
%              School of Computer Science, CMU
%    Contact : Jakramate.b@cmu.ac.th
% =========================================================================
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% =========================================================================
% INPUT:
%   funOjb : Objective function handle  (i.e. @your_function_name)
%   options: A structure of options, possible fields are
%      REQUIRED
%         options.ub
%             Upperbound on possible solution to the funObj
%         options.lb
%             Lowerbound on possible solution to the funObj
%         options.optmType
%             'min' = to minimise the funObj
%             'max' = to maximise the funObj
%      OPTIONAL
%         options.verbose
%             Print useful information during the run
%             [default = false]
%         options.maxFE
%             Maximum function evaluation allowed
%             [default = 100000]
%         options.theta %???
%             The size of best individuals retained
%             [default = 0.25 * size(x,1)] %??? I see 0.3 in the paper (and indeed, below in the code), and the meaning of theta is said to be threshold to define weakly correlated
%         options.m
%             Size of population subset
%         options.mc
%             Size of population subset for correlation analysis (WI)
%         options.c
%             Max variable group size in (SM)
%         options.eta
%             Mixing weight of new/old population
%     x:  initial population in the form of an 'n x d' matrix where 'n' is the
%         number of populations and 'd' is the dimensionality of the problem
%         (note: variable 'popSize' will be set to 'n');
%     varargin: Arguments to 'funObj'
%--------------------------------------------------------------------------
% OUTPUT:
%     optm_x : Optimal solution to the objective
%     fv     : Optimal objective value;
%     stats  : Info. about dynamics of evolution
%         stats.fv  = best fitness in each generation
%         stats.dev = std of fitnesses in the pop in each generation
%         stats.gap = distance of the best individual from x*
%     options: Options used (for follow up analysis)
% =========================================================================

% The algorithm begins at line 205.

function [optm_x fv stats options iter weak] = eda_mcc_mi(funObj, options, x, varargin)


%% preprocessing stage
[popSize d]     = size(x);
options.popSize = popSize;

% threshold for weakly dependent variables
if isfield(options, 'theta');
    theta = options.theta;
else
    theta = 0.3;  % as suggested in Dong et al.
end

if ~isfield(options, 'ub');
    error('Please provide upperbound on function');
end

if ~isfield(options, 'lb');
    error('Please provide lowerbound on function');
end

if ~isfield(options, 'optmType')
    error('Please provide optimisation type {min, max}');
end

% print some useful info?
if isfield(options, 'verbose')
    verbose = options.verbose;
else
    verbose = false;
    options.verbose = verbose;
end

% define maximum generation
if isfield(options, 'maxFE')
    maxFE = options.maxFE;
else
    maxFE = 100000;
    options.maxFE = maxFE;
end

% define sub-population size
if isfield(options, 'm')
    m = options.m;
    if m > popSize
        disp('Warning: m is larger than popSize, will use, popSize/5');
        m = 0.2 * popSize; %??? paper says 0.5*popSize
        options.m = m;
    end
else
    m = 0.2 * popSize;
    options.m = m;
end


if isfield(options, 'mc')
    mc = options.mc;
    if mc > m
        disp('Warning: mCorr is larger than m, will use, m/2');
        mc = 0.5 * m;
        options.mc = mc;
    end
else
    mc = 0.5 * m;
    options.mc = mc;
end

% define variable group size
if isfield(options, 'c')
    c = options.c;
    % fail safe check
    if c > d
        disp('Warning: group size is larger than original dimension.');
        c = d;
        options.c = c;
    end
else
    c = 3;
    options.c = c;
end


% define mixing ratio oldPop/newPop
if isfield(options, 'eta')
    eta = options.eta;
    if eta > 1
        disp('Warning: Eta is between 0 and 1');
        eta = 0.3;
        options.eta = eta;
    end
    
else
    eta = 0.3;
    options.eta = eta;
end



% all initial pop are top_x
usedFE = 0;
iter   = 1;
stats.gap=zeros(1,10000);
stats.dev=stats.gap;
stats.fv=stats.gap;


while true
 
    
    %% evaluating old pop
    remainFE  = maxFE - usedFE; 
    if remainFE < m
        break;
    end
    if remainFE == 0
        break;
    elseif  remainFE - popSize >= 0
        [fval xstar] = funObj(x, varargin{:});
        %usedFE       = usedFE + popSize;
    else
        [fval xstar] = funObj(x(1:remainFE,:),varargin{:});
        %usedFE       = usedFE + remainFE;
    end
    
    % selecting top 'topSize' individuals based on fitness values
    if strcmp(options.optmType, 'max')
        [fsort_a, a] = sort(fval, 'descend');  % maximising
    else
        [fsort_a, a] = sort(fval, 'ascend');   % minimising
    end
    optm_x = x(a(1),:);                % save optimal x
    fv     = fval(a(1));               % save optimal func. value
    
%==============================================================================%
    
    %% EDA-MCC
    new_x = zeros(popSize,d);
    
    % (1) selecting m
    %r     = randperm(popSize);
    %mSub  = x(r(1:m),:);
       mSub  = x(a(1:m),:);
    
    % (2) selecting mCorr
    rCorr = randperm(m);
    mCorr = mSub(rCorr(1:mc),:); %disp('check me...');pause
    
%     % (3) finding correlation        
%     Co    = corr(mCorr) .* ~eye(d);
%     sCo   = sum(Co <= theta);
    
     %(3) Finding mutual information
       %mCorr = (mCorr)';
%        Co = zeros(d);
%        for i = 1:d
%            for j = 1:d
%                %Co(i,j)   = MI([mCorr(:,i),mCorr(:,j)]);
%                Co(i,j)   = mutualinfo(mCorr(:,i),mCorr(:,j));
%            end
%        end
       Co = MI(mCorr);
       %Co = Co  .* ~eye(d); 
       sCo   = sum(Co <= theta);
  
    
    % finding weak correlation (i.e. less than theta)
    % and strong correlation
    weakIdx = (sCo == d);
    weak(iter) = sum(weakIdx);
    strgIdx = (sCo ~= d);
    
    % (4) WI, for weak corr fit univariate Guassian based on m
    mUni = mean(mSub(:,weakIdx));
    sUni = std(mSub(:,weakIdx));
    
    % sampling individual's WI component using mUni and sUni
    new_x(:,weakIdx) = (randn(popSize,sum(weakIdx)) ...
                        .* repmat(sUni,popSize,1)) ...        
                        +  repmat(mUni,popSize,1);
    
    % (5) SM, for strong corr fit multivariate Gaussian based on m
    sSize  = sum(strgIdx);
    nGroup = ceil(sSize/c);
    tmp    = find(strgIdx); 
    psIdx  = tmp(randperm(sSize));
    
    for i=0:nGroup-1
        lim = (i*c)+1;      
        curIdx = psIdx(lim:(min(lim+c,sSize)));
        new_x(:,curIdx) = mvnrnd(mean(mSub(:,curIdx)),...
            cov(mSub(:,curIdx)),...
            popSize);
    end
    
    %% make sure that newpop is bounded by ub and lb
    new_x(new_x>options.ub) = options.ub;
    new_x(new_x<options.lb) = options.lb;
    
    %% evaluating new pop
    remainFE  = maxFE - usedFE;
    if remainFE == 0
        break;
    elseif  remainFE - popSize >= 0
        [fval xstar] = funObj(new_x, varargin{:});
        usedFE       = usedFE + popSize;
    else
        [fval xstar] = funObj(new_x(1:remainFE,:),varargin{:});
        usedFE       = usedFE + remainFE;
    end
    
    % selecting top 'topSize' individuals based on fitness values
    if strcmp(options.optmType, 'max')
        [fsort_b, b] = sort(fval, 'descend');  % maximising
    else
        [fsort_b, b] = sort(fval, 'ascend');   % minimising
    end
    optm_nx = new_x(b(1),:);                % save optimal parameter
    fv_n    = fval(b(1));                   % save optimal func. value
    new_x(end,:) = [];
      
    %% merging new pop with the old pop
%    x = [x(a(1:ceil(length(a)*eta)),:); % at last generation length(a) can be smaller than popSize
%        new_x(b(1:floor(length(a)*(1-eta))),:)];

%>>>
% either replace the population and use elitism
    x = [x(a(1),:);... % elitism
        new_x]; 

% or merge the two populations to retain the best fitness points from both
%    all_x=[x(a(2:length(a)),:);new_x];
%    if strcmp(options.optmType, 'max')
%       [~,ab]=sort([a(2:length(a)); b],'descend'); 
%    else
%          [~,ab]=sort([a(2:length(a)); b],'ascend');
%    end	  
%     x = [x(a(1),:);... % elitism
%          all_x(ab(1:length(b)-1),:)];
%<<<	  
    
    
    if strcmp(options.optmType, 'max')
        optm_x = max(optm_nx, optm_x);  % maximising
        fv     = max(fv_n, fv);
    else
        optm_x = min(optm_nx, optm_x);  % minimising
        fv     = min(fv_n, fv);
    end
    
    
    % collecting stats
    stats.fv(iter)  = fv;
    stats.dev(iter) = std(fval);
    stats.gap(iter) = norm(xstar - optm_x);
    
    %%========================  VISUALISATION  ===========================
    if verbose
        if iter==1; figure(1); end
        % plotting function value
        plot(stats.fv(1:iter),'x','LineWidth', 2); drawnow;
        %ylim([stats.fv(iter)-10, stats.fv(iter)+10]);
        title(['Best FV:' num2str(fv)]); drawnow;
        xlabel('generations');
        ylabel('Function value');
    end
    iter = iter + 1;
end
