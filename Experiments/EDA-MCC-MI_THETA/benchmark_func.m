
function [fit optx] = benchmark_func(x, func_num, shift)
global initial_flag jrandflag jrand lb ub
persistent fhd

if (initial_flag == 0)
    % Search Range
    if (contains(func_num, [2 5 10 15]))
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
    
    switch func_num
        % Separable, D = 1000
        case 0
            fhd = str2func('sphere_shift_func');
        case 1
            fhd = str2func('elliptic_shift_func');
        case 2
            fhd = str2func('rastrigin_shift_func');
        case 3
            fhd = str2func('ackley_shift_func');
            % Single-group m-nonseparable, D = 1000, m = 50
        case 4
            fhd = str2func('elliptic_group1_shift_rot_func');
        case 5
            fhd = str2func('rastrigin_group1_shift_rot_func');
        case 6
            fhd = str2func('ackley_group1_shift_rot_func');
        case 7
            fhd = str2func('schwefel_group1_shift_func');
        case 8
            fhd = str2func('rosenbrock_group1_shift_func');
            % D/(2m)-group m-nonseparable, D = 1000, m = 50
        case 9
            fhd = str2func('elliptic_group10_shift_rot_func');
        case 10
            fhd = str2func('rastrigin_group10_shift_rot_func');
        case 11
            fhd = str2func('ackley_group10_shift_rot_func');
        case 12
            fhd = str2func('schwefel_group10_shift_func');
        case 13
            fhd = str2func('rosenbrock_group10_shift_func');
            % D/m-group m-nonseparable, D = 1000, m = 50
        case 14
            fhd = str2func('elliptic_group20_shift_rot_func');
        case 15
            fhd = str2func('rastrigin_group20_shift_rot_func');
        case 16
            fhd = str2func('ackley_group20_shift_rot_func');
        case 17
            fhd = str2func('schwefel_group20_shift_func');
        case 18
            fhd = str2func('rosenbrock_group20_shift_func');
            % Fully-nonseparable, D = 1000
        case 19
            fhd = str2func('schwefel_shift_func');
        case 20
            fhd = str2func('rosenbrock_shift_func');
	case 21 % ADDED: modified ellipse functions (3 versions)
             fhd = str2func('F3_of_cec05_a10');
	case 22 
             fhd = str2func('F3_of_cec05_a50');	     
	case 23 
             fhd = str2func('F3_of_cec05_a100');
	case 24 % ADDED: Sschwefel Function- F6 in Dong et al; unshifted version is F2 in [7] 	     
             fhd = str2func('schwefel');
    end
    
    jrandflag = 0;
    if (jrandflag == 1)
        jrand = Randomizer(func_num);
    end
end

[fit optx] = feval(fhd, x, shift);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Added F3_of_cec05 with a=10
% 	3.Shifted Rotated High Conditioned Elliptic Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx]=F3_of_cec05_a10(x,shift)
global initial_flag
persistent o M
[ps,D]=size(x);
if initial_flag==0
    load high_cond_elliptic_rot_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-100+200*rand(1,D);
    end
    c=1;
    if D==2,load elliptic_M_D2,
    elseif D==10,load elliptic_M_D10,
    elseif D==30,load elliptic_M_D30,
    elseif D==50,load elliptic_M_D50,
    else 
        A=normrnd(0,1,D,D);[M,r]=cGram_Schmidt(A);
    end
    initial_flag=1;
end
if shift==eps,
   x=x-repmat(o,ps,1); optx=o;
else
  x=x-repmat(shift,ps,D);optx=shift;
end
x=x*M;
%a=1e+6; % how about:
a=10;
fit=0;
for i=1:D
fit=fit+a.^((i-1)/(D-1)).*x(:,i).^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Added F3_of_cec05 with a=50
% 	3.Shifted Rotated High Conditioned Elliptic Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx]=F3_of_cec05_a50(x,shift)
global initial_flag
persistent o M
[ps,D]=size(x);
if initial_flag==0
    load high_cond_elliptic_rot_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-100+200*rand(1,D);
    end
    c=1;
    if D==2,load elliptic_M_D2,
    elseif D==10,load elliptic_M_D10,
    elseif D==30,load elliptic_M_D30,
    elseif D==50,load elliptic_M_D50,
    else 
        A=normrnd(0,1,D,D);[M,r]=cGram_Schmidt(A);
    end
    initial_flag=1;
end
if shift==eps,
   x=x-repmat(o,ps,1); optx=o;
else
  x=x-repmat(shift,ps,D);optx=shift;
end
x=x*M;
%a=1e+6; % how about:
a=50;
fit=0;
for i=1:D
fit=fit+a.^((i-1)/(D-1)).*x(:,i).^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Added F3_of_cec05 with a=100
% 	3.Shifted Rotated High Conditioned Elliptic Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx]=F3_of_cec05_a100(x,shift)
global initial_flag
persistent o M
[ps,D]=size(x);
if initial_flag==0
    load high_cond_elliptic_rot_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-100+200*rand(1,D);
    end
    c=1;
    if D==2,load elliptic_M_D2,
    elseif D==10,load elliptic_M_D10,
    elseif D==30,load elliptic_M_D30,
    elseif D==50,load elliptic_M_D50,
    else 
        A=normrnd(0,1,D,D);[M,r]=cGram_Schmidt(A);
    end
    initial_flag=1;
end
if shift==eps,
   x=x-repmat(o,ps,1); optx=o;
else
  x=x-repmat(shift,ps,D);optx=shift;
end
x=x*M;
%a=1e+6; % how about:
a=100;
fit=0;
for i=1:D
fit=fit+a.^((i-1)/(D-1)).*x(:,i).^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADDED: Sschwefel Function- F6 in Dong et al; unshifted version is F2 in [7] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit,optx] = schwefel(x, shift)
global initial_flag jrandflag jrand lb ub
persistent o

[ps D] = size(x);
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub);
        o = o';
        save 'datafiles/f01_o.mat' o; % there is no f00_o.mat' 
    else
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            load 'datafiles/f01_o.mat'; % there is no f24_o.mat' 
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    initial_flag = 1;
end
x    = x-repmat(o,ps,1);
%fit  = sphere_func(x);
optx = o - 1;

fit = sum( (repmat(x(:,1),1,D) - x.^2).^2 + (x-1).^2 , 2);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sphere Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = ellipsoid_func(x, A)

fit = sum(x*A*x,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sphere Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = sphere_func(x)

fit = sum(x.*x, 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elliptic Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = elliptic_func(x)

[ps, D] = size(x);
a = 1e+6;
fit = 0;
for i = 1:D
    fit = fit+a.^((i-1)/(D-1)).*x(:,i).^2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotated Elliptic Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = elliptic_rot_func(x, M)

x = x*M;
fit = elliptic_func(x);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schwefel's Problem 1.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = schwefel_func(x)

[ps D] = size(x);
fit = 0;
for i = 1:D
    fit = fit + sum(x(:,1:i),2).^2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rosenbrock's Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = rosenbrock_func(x)

[ps D] = size(x);
fit = sum(100.*(x(:,1:D-1).^2-x(:,2:D)).^2+(x(:,1:D-1)-1).^2,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rastrigin's Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = rastrigin_func(x)

fit = sum(x.*x-10*cos(2*pi*x)+10, 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotated Rastrigin's Fucntion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = rastrigin_rot_func(x, M)

x = x*M;
fit = rastrigin_func(x);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ackley's Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = ackley_func(x)

[ps, D] = size(x);
fit = sum(x.^2,2);
fit = 20-20.*exp(-0.2.*sqrt(fit./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotated Ackley's Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = ackley_rot_func(x, M)

x = x*M;
fit = ackley_func(x);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F0: Shifted Sphere Function
% D = 1000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = sphere_shift_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o

[ps D] = size(x);
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub);
        o = o';
        save 'datafiles/f01_o.mat' o; % there is no f00_o.mat' 
    else
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            load 'datafiles/f01_o.mat'; % there is no f00_o.mat' 
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    initial_flag = 1;
end
x    = x-repmat(o,ps,1);
fit  = sphere_func(x);
optx = o;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F1: Shifted Elliptic Function
% D = 1000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = elliptic_shift_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o

[ps D] = size(x);
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub);
        o = o';
        save 'datafiles/f01_o.mat' o;
    else
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            load 'datafiles/f01_o.mat';
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    initial_flag = 1;
end

x    = x-repmat(o,ps,1);
fit  = elliptic_func(x);
optx = o;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F2: Shifted Rastrigin's Function
% D = 1000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = rastrigin_shift_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o

[ps D] = size(x);
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub);
        o = o';
        save 'datafiles/f02_o.mat' o;
    else
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            load 'datafiles/f02_o.mat';
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    initial_flag = 1;
end

x = x-repmat(o,ps,1);
fit = rastrigin_func(x);
optx = o;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F3: Shifted Ackley's Function
% D = 1000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = ackley_shift_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o

[ps D] = size(x);
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub);
        o = o';
        save 'datafiles/f03_o.mat' o;
    else
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            load 'datafiles/f03_o.mat';
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    initial_flag = 1;
end

x = x-repmat(o,ps,1);
fit = ackley_func(x);
optx = o;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F4: Single Group Shifted and Rotated Elliptic Function
% D = 1000, m = 50
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = elliptic_group1_shift_rot_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o p M

[ps D] = size(x);
m = 50;

if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub);
        o = o';
        p = jrand.createPermVector(D);
        p = p'+1;
        M = jrand.createRotMatrix(m);
        save 'datafiles/f04_opm.mat' o p M;
    else
        load 'datafiles/f04_opm.mat';
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    if (D ~= 1000)
        %disp('F4 error: only support D = 1000 now');
        % exit(4);
        idx = p<D;
        p   = p(idx);
    end
    initial_flag = 1;
end

a = 1e+6;
x = x-repmat(o,ps,1);
fit = a*elliptic_rot_func(x(:,p(1:m)), M) + elliptic_func(x(:,p((m+1):end)));
optx = o;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F5: Single Group Shifted and Rotated Rastrigin's Function
% D = 1000, m = 50
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = rastrigin_group1_shift_rot_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o p M

[ps D] = size(x);
m = 50;
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub);
        o = o';
        p = jrand.createPermVector(D);
        p = p'+1;
        M = jrand.createRotMatrix(m);
        save 'datafiles/f05_opm.mat' o p M;
    else
        load 'datafiles/f05_opm.mat';
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    if (D ~= 1000)
        %disp('F5 error: only support D = 1000 now');
        %exit(5);
        idx = p<D;
        p   = p(idx);
    end
    initial_flag = 1;
end

a = 1e+6;
x = x-repmat(o,ps,1);
fit = a*rastrigin_rot_func(x(:,p(1:m)), M) + rastrigin_func(x(:,p((m+1):end)));
optx = o;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F6: Single Group Shifted and Rotated Ackley's Function
% D = 1000, m = 50
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = ackley_group1_shift_rot_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o p M

[ps D] = size(x);
m = 50;
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub);
        o = o';
        p = jrand.createPermVector(D);
        p = p'+1;
        M = jrand.createRotMatrix(m);
        save 'datafiles/f06_opm.mat' o p M;
    else
        load 'datafiles/f06_opm.mat';
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    if (D ~= 1000)
        %disp('F6 error: only support D = 1000 now');
        %exit(6);
        idx = p<D;
        p   = p(idx);
    end
    initial_flag = 1;
end

a = 1e+6;
x = x-repmat(o,ps,1);
fit = a*ackley_rot_func(x(:,p(1:m)), M) + ackley_func(x(:,p((m+1):end)));
optx = o;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F7: Single Group Shifted Schwefel's Problem 1.2
% D = 1000, m = 50
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = schwefel_group1_shift_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o p

[ps D] = size(x);
m = 50;
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub);
        o = o';
        p = jrand.createPermVector(D);
        p = p'+1;
        save 'datafiles/f07_op.mat' o p;
    else
        load 'datafiles/f07_op.mat';
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    if (D ~= 1000)
        %disp('F7 error: only support D = 1000 now');
        %exit(7);
        idx = p<D;
        p   = p(idx);
    end
    initial_flag = 1;
end

a = 1e+6;
x = x-repmat(o,ps,1);
fit = a*schwefel_func(x(:,p(1:m))) + sphere_func(x(:,p((m+1):end)));
optx = o;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F8: Single Group Shifted Rosenbrock's Function
% D = 1000, m = 50
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = rosenbrock_group1_shift_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o p

[ps D] = size(x);
m = 50;
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub-1);
        o = o';
        p = jrand.createPermVector(D);
        p = p'+1;
        save 'datafiles/f08_op.mat' o p;
    else
        load 'datafiles/f08_op.mat';
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    if (D ~= 1000)
        %disp('F8 error: only support D = 1000 now');
        %exit(8);
        idx = p<D;
        p   = p(idx);
    end
    initial_flag = 1;
end

a = 1e+6;
x = x-repmat(o,ps,1);
fit = a*rosenbrock_func(x(:,p(1:m))) + sphere_func(x(:,p((m+1):end)));
optx = o;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F9: D/(2m)-group Shifted and Rotated Elliptic Function
% D = 1000, m = 50, D/(2m) = 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = elliptic_group10_shift_rot_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o p M

[ps D] = size(x);
m = 50;
G = D/m/2;
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub);
        o = o';
        p = jrand.createPermVector(D);
        p = p'+1;
        M = jrand.createRotMatrix(m);
        save 'datafiles/f09_opm.mat' o p M;
    else
        load 'datafiles/f09_opm.mat';
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    if (D ~= 1000)
        %disp('F9 error: only support D = 1000 now');
        %exit(9);
        idx = p<D;
        p   = p(idx);
    end
    initial_flag = 1;
end

x = x-repmat(o,ps,1);
fit = 0;
for k = 1:G
    index = ((k-1)*m+1):(k*m);
    fit = fit + elliptic_rot_func(x(:,p(index)), M);
end
fit = fit + elliptic_func(x(:,p((G*m+1):end)));
optx = o;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F10: D/(2m)-group Shifted and Rotated Rastrigin's Function
% D = 1000, m = 50, D/(2m) = 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = rastrigin_group10_shift_rot_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o p M

[ps D] = size(x);
m = 50;
G = D/m/2;
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub);
        o = o';
        p = jrand.createPermVector(D);
        p = p'+1;
        M = jrand.createRotMatrix(m);
        save 'datafiles/f10_opm.mat' o p M;
    else
        load 'datafiles/f10_opm.mat';
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    if (D ~= 1000)
        %disp('F10 error: only support D = 1000 now');
        %exit(10);
        idx = p<D;
        p   = p(idx);
    end
    initial_flag = 1;
end

x = x-repmat(o,ps,1);
fit = 0;
for k = 1:G
    index = ((k-1)*m+1):(k*m);
    fit = fit + rastrigin_rot_func(x(:,p(index)), M);
end
fit = fit + rastrigin_func(x(:,p((G*m+1):end)));
optx = o;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F11: D/(2m)-group Shifted and Rotated Ackley's Function
% D = 1000, m = 50, D/(2m) = 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = ackley_group10_shift_rot_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o p M

[ps D] = size(x);
m = 50;
G = D/m/2;
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub);
        o = o';
        p = jrand.createPermVector(D);
        p = p'+1;
        M = jrand.createRotMatrix(m);
        save 'datafiles/f11_opm.mat' o p M;
    else
        load 'datafiles/f11_opm.mat';
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    if (D ~= 1000)
        %disp('F11 error: only support D = 1000 now');
        %exit(11);
        idx = p<D;
        p   = p(idx);
    end
    initial_flag = 1;
end

x = x-repmat(o,ps,1);
fit = 0;
for k = 1:G
    index = ((k-1)*m+1):(k*m);
    fit = fit + ackley_rot_func(x(:,p(index)), M);
end
fit = fit + ackley_func(x(:,p((G*m+1):end)));
optx = o;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F12: D/(2m)-group Shifted Schwefel's Problem 1.2
% D = 1000, m = 50, D/(2m) = 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = schwefel_group10_shift_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o p

[ps D] = size(x);
m = 50;
G = D/m/2;
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub);
        o = o';
        p = jrand.createPermVector(D);
        p = p'+1;
        save 'datafiles/f12_op.mat' o p;
    else
        load 'datafiles/f12_op.mat';
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    if (D ~= 1000)
        %disp('F12 error: only support D = 1000 now');
        %exit(12);
        idx = p<D;
        p   = p(idx);
    end
    initial_flag = 1;
end

x = x-repmat(o,ps,1);
fit = 0;
for k = 1:G
    index = ((k-1)*m+1):(k*m);
    fit = fit + schwefel_func(x(:,p(index)));
end
fit = fit + sphere_func(x(:,p((G*m+1):end)));
optx = o;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F13: D/(2m)-group Shifted Rosenbrock's Function
% D = 1000, m = 50, D/(2m) = 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = rosenbrock_group10_shift_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o p

[ps D] = size(x);
m = 50;
G = D/m/2;
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub-1);
        o = o';
        p = jrand.createPermVector(D);
        p = p'+1;
        save 'datafiles/f13_op.mat' o p;
    else
        load 'datafiles/f13_op.mat';
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    if (D ~= 1000)
        %disp('F13 error: only support D = 1000 now');
        %exit(13);
        idx = p<D;
        p   = p(idx);
    end
    initial_flag = 1;
end

x = x-repmat(o,ps,1);
fit = 0;
for k = 1:G
    index = ((k-1)*m+1):(k*m);
    fit = fit + rosenbrock_func(x(:,p(index)));
end
fit = fit + sphere_func(x(:,p((G*m+1):end)));
optx = o;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F14: D/m-group Shifted and Rotated Elliptic Function
% D = 1000, m = 50, D/m = 20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = elliptic_group20_shift_rot_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o p M

[ps D] = size(x);
m = 50;
G = D/m;
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub);
        o = o';
        p = jrand.createPermVector(D);
        p = p'+1;
        M = jrand.createRotMatrix(m);
        save 'datafiles/f14_opm.mat' o p M;
    else
        load 'datafiles/f14_opm.mat';
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    if (D ~= 1000)
        %disp('F14 error: only support D = 1000 now');
        %exit(14);
        idx = p<D;
        p   = p(idx);
    end
    initial_flag = 1;
end

x = x-repmat(o,ps,1);
fit = 0;
for k = 1:G
    index = ((k-1)*m+1):(k*m);
    fit = fit + elliptic_rot_func(x(:,p(index)), M);
end
optx = o;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F15: D/m-group Shifted and Rotated Rastrigin's Function
% D = 1000, m = 50, D/m = 20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = rastrigin_group20_shift_rot_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o p M

[ps D] = size(x);
m = 50;
G = D/m;
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub);
        o = o';
        p = jrand.createPermVector(D);
        p = p'+1;
        M = jrand.createRotMatrix(m);
        save 'datafiles/f15_opm.mat' o p M;
    else
        load 'datafiles/f15_opm.mat';
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    if (D ~= 1000)
        %disp('F15 error: only support D = 1000 now');
        %exit(15);
        idx = p<D;
        p   = p(idx);
    end
    initial_flag = 1;
end

x = x-repmat(o,ps,1);
fit = 0;
for k = 1:G
    index = ((k-1)*m+1):(k*m);
    fit = fit + rastrigin_rot_func(x(:,p(index)), M);
end
optx = o;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F16: D/m-group Shifted and Rotated Ackley's Function
% D = 1000, m = 50, D/m = 20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = ackley_group20_shift_rot_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o p M

[ps D] = size(x);
m = 50;
G = D/m;
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub);
        o = o';
        p = jrand.createPermVector(D);
        p = p'+1;
        M = jrand.createRotMatrix(m);
        save 'datafiles/f16_opm.mat' o p M;
    else
        load 'datafiles/f16_opm.mat';
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    if (D ~= 1000)
        %disp('F16 error: only support D = 1000 now');
        %exit(16);
        idx = p<D;
        p   = p(idx);
    end
    initial_flag = 1;
end

x = x-repmat(o,ps,1);
fit = 0;
for k = 1:G
    index = ((k-1)*m+1):(k*m);
    fit = fit + ackley_rot_func(x(:,p(index)), M);
end
optx = o;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F17: D/m-group Shifted Schwefel's Problem 1.2
% D = 1000, m = 50, D/m = 20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = schwefel_group20_shift_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o p

[ps D] = size(x);
m = 50;
G = D/m;
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub);
        o = o';
        p = jrand.createPermVector(D);
        p = p'+1;
        save 'datafiles/f17_op.mat' o p;
    else
        load 'datafiles/f17_op.mat';
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    if (D ~= 1000)
        %disp('F17 error: only support D = 1000 now');
        %exit(17);
        idx = p<D;
        p   = p(idx);
    end
    initial_flag = 1;
end

x = x-repmat(o,ps,1);
fit = 0;
for k = 1:G
    index = ((k-1)*m+1):(k*m);
    fit = fit + schwefel_func(x(:,p(index)));
end
optx = o;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F18: D/m-group Shifted Rosenbrock's Function
% D = 1000, m = 50, D/m = 20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = rosenbrock_group20_shift_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o p

[ps D] = size(x);
m = 50;
G = D/m;
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub-1);
        o = o';
        p = jrand.createPermVector(D);
        p = p'+1;
        save 'datafiles/f18_op.mat' o p;
    else
        load 'datafiles/f18_op.mat';
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    if (D ~= 1000)
        %disp('F18 error: only support D = 1000 now');
        %exit(18);
        idx = p<D;
        p   = p(idx);
    end
    initial_flag = 1;
end

x = x-repmat(o,ps,1);
fit = 0;
for k = 1:G
    index = ((k-1)*m+1):(k*m);
    fit = fit + rosenbrock_func(x(:,p(index)));
end
optx = o;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F19: Shifted Schwefel's Problem 1.2
% D = 1000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = schwefel_shift_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o

[ps D] = size(x);
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub);
        o = o';
        save 'datafiles/f19_o.mat' o;
    else
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            load 'datafiles/f19_o.mat';
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    if (D ~= 1000)
        %disp('F19 error: only support D = 1000 now');
        %exit(19);
    end
    initial_flag = 1;
end
optx = o;

x = x-repmat(o,ps,1);
fit = schwefel_func(x);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F20: Shifted Rosenbrock's Function
% D = 1000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = rosenbrock_shift_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o

[ps D] = size(x);
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub-1);
        o = o';
        save 'datafiles/f20_o.mat' o;
    else
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            load 'datafiles/f20_o.mat';
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    if (D ~= 1000)
        %disp('F20 error: only support D = 1000 now');
        %exit(20);
    end
    initial_flag = 1;
end

x = x-repmat(o,ps,1);
fit = rosenbrock_func(x);
optx = o;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F21: Shifted Ellipsoid Function
% D = 1000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit optx] = ellipsoid_shift_func(x,shift)
global initial_flag jrandflag jrand lb ub
persistent o A

[ps D] = size(x);
if (initial_flag == 0)
    if (jrandflag == 1)
        o = jrand.createShiftVector(D, lb, ub);
        o = o';
        save 'datafiles/f21_o.mat' o;
    else
        % eps shift means using theirs
        % predefined shifts
        if shift == eps
            load 'datafiles/f21_o.mat';
            o = o(1:D);
        else
            o = zeros(1,D) + shift;
        end
    end
    initial_flag = 1;
end
x    = x-repmat(o,ps,1);
fit  = ellipsoid_func(x, A);
optx = o;

% classical Gram Schmid 
 function [q,r] = cGram_Schmidt (A)
% computes the QR factorization of $A$ via
% classical Gram Schmid 
% 
 [n,m] = size(A); 
 q = A;    
 for j=1:m
     for i=1:j-1 
         r(i,j) = q(:,j)'*q(:,i);
     end
     for i=1:j-1   
       q(:,j) = q(:,j) -  r(i,j)*q(:,i);
     end
     t =  norm(q(:,j),2 ) ;
     q(:,j) = q(:,j) / t ;
     r(j,j) = t  ;
 end 
