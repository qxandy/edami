function run_edaMCC_MI(pD,pNp,pFun,pRuns)
%addpath AMaLGaM
addpath CEC-Benchmark
%addpath CMA-ES      
addpath EDA-MCC
%addpath RP-EDA
%addpath Vanilla-EDA
%addpath RP-EDA/libs/mtimesx_20110223
%addpath RP-EDA/libs/
addpath CEC-Benchmark/datafiles

% for experimental purpose
%s = RandStream('mcg16807','Seed',121121); RandStream.setDefaultStream(s);

D    = pD;         % dimensionality of the problem
NP   = pNp;%300;          % population size
runs = pRuns; %25;            % repetitions
global initial_flag; % the global flag used in test suite
FE= 1*10^4*D; %3*10^6; % budget of function evaluations
FUNCTIONS=pFun; %0 2 3 5 6 8 10 11 13 15 16 18 20]; % 0=sphere, the rest are the multimodal functions from CEC10 benchmark
METHODS={'mcc'}; %METHOD %'rp','cma','mcc','amalgam'
Verbose=0;

for func_num = FUNCTIONS
    
    if (contains(func_num, [2 5 10 15])) % if the func_num is one of the numbers in the array -> true
        lb = -5;
        ub = 5;
    else
        if (contains(func_num, [3 6 11 16]))
            lb = -32;
            ub = 32;
        else
            lb = -100;
            ub = 100;
        end
    end
    
    for run = 1:runs
        close
        % should set the flag to 0 for each run, each function
        initial_flag = 0;
        pop = lb + (ub-lb) * rand(NP, D); % initial population of NP individuals
                       
        % ---- RP_EDA ----------------------------------------------------------
todo=0; for mth=1:length(METHODS), if strcmp(METHODS{mth},'rp'),todo=1;end;end
if todo, disp('Running RP-EDA...')
        opt1.maxFE    = FE; 
        opt1.optmType = 'min';
        opt1.ub       = ub;
        opt1.lb       = lb;
        opt1.topSize  = round(NP/4);       
        opt1.k        = 3;
        opt1.rpmSize  = D;
        opt1.verbose  = Verbose;
        
        t1 = cputime;
        % set the last argument (shift amount) to 'eps' to use the suite's
        % predefined shift values otherwise the function will be shifted by
        % 'shift amount'
        [x1, val1(run), stats1, opts1] = rp_eda(@benchmark_func, opt1, pop, func_num, eps);
        % print & save
        fprintf(1, 'RP-EDA Elapsed CPU time: %f\n', cputime-t1);
        fprintf(1, 'func_num = %d, run = %d\n', func_num, run);
        fprintf(1, 'min(val) RP= %f\n\n', val1(run));
	Time1=cputime-t1;
        save(['results-ours-3e06FE/rp_F' num2str(func_num) '_rep_' num2str(run)],'x1','val1','stats1','opts1','Time1');   
end        
        
        %% ---- EDA-MCC --------------------- this is the stuff to run ------------------------------------
todo=0; for mth=1:length(METHODS), if strcmp(METHODS{mth},'mcc'),todo=1;end;end
if todo, disp('Running EDA-MCC-MI...')
        opt2.maxFE    = FE;
        opt2.optmType = 'min';
        opt2.ub       = ub;
        opt2.lb       = lb;
        opt2.theta    = 0.3;    % weakly dependent threshold (theta in the paper)
        opt2.m        = round(NP/2); % (m in the paper) nr of selected individuals, where NP=population size
        opt2.mc       = 100;    % mCorr size (m_corr in the paper)
        opt2.c        = min(round(D/5),round(NP/15)); %max group size (c in the paper). Note a typo on page 811 of the MCC paper, for their 500-dimensional experiments 
	  % where they say they have Nr selected points = 100, block size c = 100 must be wrong and doesn't work.
          % Indeed this setting would mean to estimate a 100 x 100 covariance from 100 points -- bad eigenvalue misestimation would occur.
          % We either need to go lower with the block size (as here, c=20) or else one needs to increase the population size NP 
	  % With these settings we have c x c = 20 x 20 covariance blocks to estimate from opt2.m = round(NP/2) points
        opt2.eta      = 0;    % how many old pop are kept [0 1] % not used any more in eda_mcc.m
        opt2.verbose  = Verbose; % 0 or 1
        tic;
        [x2, val2(run), stats2, opts2 iter] = eda_mcc_mi(@benchmark_func, opt2, pop, func_num, eps); % eps means to use the default shifts in the CEC2010 benchmark functions
	    % otherwise you can specify a shift (i.e. a number) for the last parameter 
        % print & save
        ti=toc;
        fprintf(1, 'EDA-MCC-MI Elapsed CPU time: %f\n', ti); 
        fprintf(1, 'EDA-MCC-MI iterations: %d\n', iter); 
        fprintf(1, 'func_num = %d, run = %d\n', func_num, run);
        fprintf(1, 'min(val) MCC=%f\n\n',val2(run));            
	%Time2=cputime-t2;
    	B(run,:) = stats2;
        save(['BFG_'  num2str(D) '_' num2str(NP) '_F' num2str(func_num)],'B');
        %FV(run,:) = stats2.fv;
        %save(['BFV_'  num2str(D) '_' num2str(NP) '_F' num2str(func_num)],'FV');
        
        %save(['mcc_F' num2str(func_num) '_rep_' num2str(run)],'x2','val2','stats2','opts2','Time2');   
end

%=== END ================== END ===================================================================================%

        %% ---- CMA-ES ---------------------------------------------------------
todo=0; for mth=1:length(METHODS), if strcmp(METHODS{mth},'cma'),todo=1;end;end
if todo, disp('Running CMA-ES...')
        initial_flag = 0;        
        opt3.MaxFunEvals   = FE;
	%opt3.StopFitness=-1;
        opt3.LBounds       = lb;
        opt3.UBounds       = ub;
        opt3.SaveVariables = 'off';
        opt3.DispFinal     = 'off';
        opt3.LogTime       = 0;
        opt3.DispModulo    = Inf;       % print at every ... iterations Inf = no print
	%opt3.DispModulo    = 1;
        t3 = cputime;
        [x3, fv, counteval, stopflag, stats3, bestever,trajectory] ... %SHOULD WE PLUG pop OR CENTRE OF SEARCH SPACE?
        = cmaes('benchmark_func_slim', pop', (ub-lb)/3, opt3, func_num, eps);
        % print & save
	val3(run)=bestever.f;stats3.fv=trajectory;
        fprintf(1, 'CMA-ES Elapsed CPU time: %f\n', cputime-t3);        
        fprintf(1, 'func_num = %d, run = %d\n', func_num, run);
        fprintf(1, 'min(val) CMA-ES=%f\n',val3(run));            
	Time3=cputime-t3;
        save(['results/cma_F' num2str(func_num) '_rep_' num2str(run)],'x3','val3','stats3','Time3','counteval','stopflag','bestever');
end
        
        %% --- AMaLGaM ----------------------------------------------------------
todo=0; for mth=1:length(METHODS), if strcmp(METHODS{mth},'amalgam'),todo=1;end;end
if todo, disp('Running Amalgam...')
        opt4.LBounds = lb;
        opt4.UBounds = ub;
        opt4.rotation = 0;
        % tau pop dmd srt imp will be set by -g flag
        opt4.tau = 0.1;  
        opt4.pop = NP/3;
        opt4.dmd = 0.1;
        opt4.srt = 0.01;
        opt4.imp = 20;
        opt4.nop = 5;
        opt4.eva = FE;
        opt4.vtr = 0;%0.0000000001;
        opt4.tol = 0; %0.00001;
        
        if func_num == 20
            func_n = 7;
        else
            func_n = func_num;
        end
        
        str = sprintf('./AMaLGaM/AMaLGaM-Full -g %d %d %d %d %f %f %d %d %f %f %d %f %f %f',...
            func_n, D, opt4.LBounds, opt4.UBounds, opt4.rotation, ...
            opt4.tau, opt4.pop, opt4.nop, opt4.dmd, opt4.srt, opt4.eva,...
            opt4.vtr, opt4.imp, opt4.tol);  
        t4 = cputime;
        system(str);
        fprintf(1, 'AMaLGaM Elapsed CPU time: %f\n', cputime-t4);    
        amfull = load('best_generation_final.dat');
        x4 = amfull(1:end-2); 
        val4(run) = amfull(end-1);
        fprintf(1, 'func_num = %d, run = %d\n', func_num, run);
        fprintf(1, 'min(val) AMaLGaM=%f\n', val4(run));                    
	Time4=cputime-t4;
        save(['amalgam_F' num2str(func_num) '_rep_' num2str(run)],'x4','val4','Time4','opt4');
end

    end
end
