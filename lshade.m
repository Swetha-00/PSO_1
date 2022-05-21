%%%%%%%%%%%%%%%%%%%
%% This package is a MATLAB/Octave source code of L-SHADE which is an improved version of SHADE 1.1.
%% Note that this source code is transferred from the C++ source code version.
%% About L-SHADE, please see following papers:
%%
%% * Ryoji Tanabe and Alex Fukunaga: Improving the Search Performance of SHADE Using Linear Population Size Reduction,  Proc. IEEE Congress on Evolutionary Computation (CEC-2014), Beijing, July, 2014.
%%
%% For this package, we downloaded JADE's source code from Dr. Q. Zhang's website (http://dces.essex.ac.uk/staff/qzhang) and modified it.
%%
%% Update
%% 9/Oct/2014: incompatibilities-fix between Octave and Matlab and bug-fix of population size reduction procedure (thanks to Dr. Elsayed)
%%%%%%%%%%%%%%%%%%% 

format long;
format compact;
global wind_farm farm_power Efficiency cost_per_kW;

    pop_size = 200;
    problem_size = 100;
    Var_max = 1; % Max variant value
    Var_min = 0; % Min variant value
    X_max = ones(1,problem_size)*Var_max;
    X_min = ones(1,problem_size)*Var_min;
    X = [X_min ; X_max];
    U = zeros(1,problem_size);
    Y = zeros(pop_size,1);
    conv = zeros(2000,3); %store convergence of matrix
    fitness = zeros(); 
    max_nfes = 300 * problem_size;
tic
rand('seed', sum(100 * clock)); %create 15 digit after decimal random number

val_2_reach = 10^(-8);

Func=@Analyse_Grid;
% PSO Parameters
w=1;            % Inertia Weight
wdamp=0.99;     % Inertia Weight Damping Ratio
c1=1.5;         % Personal Learning Coefficient
c2=2.0;         % Global Learning Coefficient


for func = 1 : 1
 % optimum = 0;

  %% Record the best results
  %outcome = zeros(); 

  fprintf('\n-------------------------------------------------------\n')
  fprintf('Function = %d, Dimension size = %d\n', func, problem_size) 
  for run_id = 1 : 1
    %%  parameter settings for L-SHADE
    p_best_rate = 0.11;
    arc_rate = 1.4;
    memory_size = 5;
%     pop_size = 18 * problem_size;

    max_pop_size = pop_size;
    min_pop_size = 4.0;

    %% Initialize the main population
    % create population of [matrix pop_size x problem_size] with values
 
    popold = repmat(X(1, :),pop_size,1)+rand(pop_size, problem_size).*(repmat(X(2, :)- X(1, :),pop_size,1));
    %popold = repmat(X(1, :),pop_size,1)+randsrc(pop_size, problem_size,[0 1/3 2/3 1;0.5 0.17 0.17 0.16]).*(repmat(X(2, :)- X(1, :),pop_size,1));
    %popold = repmat(X(1, :),pop_size,1)+randsrc(pop_size, problem_size,[0 0.5 1;0.5 0.25 0.25]).*(repmat(X(2, :)- X(1, :),pop_size,1));
    
    pop = popold; % the old population becomes the current population

    for jj = 1: pop_size
        fitness(jj) = Func(pop(jj,:));
    end
    fitness = fitness';

    nfes = 0;
    bsf_fit_var = 1e+30;
    bsf_solution = zeros(1, problem_size);

    %%%%%%%%%%%%%%%%%%%%%%%% for out
    for i = 1 : pop_size
      nfes = nfes + 1;

      if fitness(i) < bsf_fit_var
	bsf_fit_var = fitness(i);
	bsf_solution = pop(i, :);
      end

	  %% if mod(nfes, 1000) == 0	
	  %% bsf_error_var = bsf_fit_var  - optimum;
	  %% if bsf_error_var < val_2_reach; bsf_error_var = 0; end;
	  %%	fprintf(sprintf('%1.16e \n', bsf_error_var));
	  %%    fprintf(sprintf('%d %1.16e \n', nfes, bsf_error_var));
	  %% end
	  
%%      if nfes > max_nfes; exit(1); end
      if nfes > max_nfes; break; end
    end
    %%%%%%%%%%%%%%%%%%%%%%%% for out

    memory_sf = 0.5 .* ones(memory_size, 1);
    memory_cr = 0.5 .* ones(memory_size, 1);
    memory_pos = 1;

    archive.NP = arc_rate * pop_size; % the maximum size of the archive
    archive.pop = zeros(0, problem_size); % the solutions stored in te archive
    archive.funvalues = zeros(0, 1); % the function value of the archived solutions

%PSO Parameter Setting Start    
 vel=zeros(pop_size,problem_size);
 pbest_pos=zeros(pop_size,problem_size);
 pbest_val = zeros(pop_size,1);
 GlobalBest_Cost=inf;
  for i=1:pop_size
     pbest_pos(i,:)=pop(i,:);
     pbest_val(i,1)=fitness(i);
     
      % Update Global Best
    if pbest_val(i,1)<GlobalBest_Cost
        
        bsf_solution=pbest_pos(i,:);
        GlobalBest_Cost=pbest_val(i,1);
    end
  end
 %END PSO Parameter Setting
    %% main loop
    index = 1;
    while nfes < max_nfes
        %*****************************PSO STARTING************************
    newpop=pop;
    for i=1:pop_size
        % Update Velocity
        vel(i,:) = w*vel(i,:)...
            +c1*rand(1,problem_size).*(pbest_pos(i,:)-pop(i,:)) ...
            +c2*rand(1,problem_size).*(bsf_solution-pop(i,:));
        % Update Position
        newpop(i,:) = newpop(i,:) + vel(i,:);
    
    
    newpop = boundConstraint(newpop,pop,X_max,X_min);
    newval = zeros(pop_size,1);            
   
      [newval(i,1)] = Func(newpop(i,:));
          
   nfes = nfes + 1;
   
    
        % Update Personal Best
        if newval(i,1)<pbest_val(i,1)
            
            pbest_pos(i,:)=newpop(i,:);
            pbest_val(i,1)=newval(i,1);
            
            % Update Global Best
            if pbest_val(i,1)<GlobalBest_Cost
                
                bsf_solution=pbest_pos(i,:);
                
            end
            
        end    
    end
    
    w=w*wdamp;
    % nfeval=nfeval+(pop_size);
     %*****************************PSO ENDING******************
     
     
%       pop = popold; % the old population becomes the current population
%       [temp_fit, sorted_index] = sort(fitness, 'ascend');
% 
%       mem_rand_index = ceil(memory_size * rand(pop_size, 1));
%       mu_sf = memory_sf(mem_rand_index);
%       mu_cr = memory_cr(mem_rand_index);
% 
%       %% for generating crossover rate
%       cr = normrnd(mu_cr, 0.1);
%       term_pos = find(mu_cr == -1);
%       cr(term_pos) = 0;
%       cr = min(cr, 1);
%       cr = max(cr, 0);
% 
%       %% for generating scaling factor
%       sf = mu_sf + 0.1 * tan(pi * (rand(pop_size, 1) - 0.5));
%       pos = find(sf <= 0);
% 
%       while ~ isempty(pos)
% 	sf(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
% 	pos = find(sf <= 0);
%       end
% 
%       sf = min(sf, 1); 
%       
%       r0 = (1 : pop_size);
%       popAll = [pop; archive.pop];
%       [r1, r2] = gnR1R2(pop_size, size(popAll, 1), r0);
%       
%       pNP = max(round(p_best_rate * pop_size), 2); %% choose at least two best solutions
%       randindex = ceil(rand(1, pop_size) .* pNP); %% select from [1, 2, 3, ..., pNP]
%       randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
%       pbest = pop(sorted_index(randindex), :); %% randomly choose one of the top 100p% solutions
% 
%       vi = pop + sf(:, ones(1, problem_size)) .* (pbest - pop + pop(r1, :) - popAll(r2, :));
%       vi = boundConstraint(vi, pop, X_max, X_min);
% 
%       mask = rand(pop_size, problem_size) > cr(:, ones(1, problem_size)); % mask is used to indicate which elements of ui comes from the parent
%       rows = (1 : pop_size)'; cols = floor(rand(pop_size, 1) * problem_size)+1; % choose one position where the element of ui doesn't come from the parent
%       jrand = sub2ind([pop_size problem_size], rows, cols); mask(jrand) = false;
%       ui = vi; ui(mask) = pop(mask);
% 
% %      children_fitness = feval(fhd, ui', func);
%       children_fitness = zeros();
%       for jj = 1: pop_size
%       children_fitness(jj) = Func(ui(jj,:));
%       end 
%       children_fitness = children_fitness';

      %%%%%%%%%%%%%%%%%%%%%%%% for out
%       for i = 1 : pop_size
%         nfes = nfes + 1;
% 
% 	if children_fitness(i) < bsf_fit_var
% 	  bsf_fit_var = children_fitness(i);
% 	  bsf_solution = ui(i, :);
% 	end

	  %% if mod(nfes, 1000) == 0	
	  %% bsf_error_var = bsf_fit_var  - optimum;
	  %% if bsf_error_var < val_2_reach; bsf_error_var = 0; end;
	  %%       fprintf(sprintf('%1.16e \n', bsf_error_var));
	  %%       fprintf(sprintf('%d %1.16e \n', nfes, bsf_error_var));
	  %%end

	%%	if nfes > max_nfes; exit(1); end
	if nfes > max_nfes; break; end
      end
      %%%%%%%%%%%%%%%%%%%%%%%% for out

%       dif = abs(fitness - children_fitness);


      %% I == 1: the parent is better; I == 2: the offspring is better
%       I = (fitness > children_fitness);
%       goodCR = cr(I == 1);  
%       goodF = sf(I == 1);
%       dif_val = dif(I == 1);
% 
% %      isempty(popold(I == 1, :))   
%       archive = updateArchive(archive, popold(I == 1, :), fitness(I == 1));
% 
%       [fitness, I] = min([fitness, children_fitness], [], 2);
%       
%       popold = pop;
%       popold(I == 2, :) = ui(I == 2, :);
% 
%       num_success_params = numel(goodCR);
% 
%       if num_success_params > 0 
% 	sum_dif = sum(dif_val);
% 	dif_val = dif_val / sum_dif;

	%% for updating the memory of scaling factor 
% 	memory_sf(memory_pos) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);
% 
% 	%% for updating the memory of crossover rate
% 	if max(goodCR) == 0 || memory_cr(memory_pos)  == -1
% 	  memory_cr(memory_pos)  = -1;
% 	else
% 	  memory_cr(memory_pos) = (dif_val' * (goodCR .^ 2)) / (dif_val' * goodCR);
% 	end
% 
% 	memory_pos = memory_pos + 1;
% 	if memory_pos > memory_size;  memory_pos = 1; end
%       end
% 
%       %% for resizing the population size
%       plan_pop_size = round((((min_pop_size - max_pop_size) / max_nfes) * nfes) + max_pop_size);
% 
%     if pop_size > plan_pop_size
% 	reduction_ind_num = pop_size - plan_pop_size;
% 	if pop_size - reduction_ind_num <  min_pop_size; reduction_ind_num = pop_size - min_pop_size;end
% 
% 	pop_size = pop_size - reduction_ind_num;
% 	for r = 1 : reduction_ind_num
% 	  [valBest,indBest] = sort(fitness, 'ascend');
% 	  worst_ind = indBest(end);
% 	  popold(worst_ind,:) = [];
% 	  pop(worst_ind,:) = [];
% 	  fitness(worst_ind,:) = [];
% 	end
% 	  
% 	archive.NP = round(arc_rate * pop_size); 
% 
% 	if size(archive.pop, 1) > archive.NP 
% 	  rndpos = randperm(size(archive.pop, 1));
% 	  rndpos = rndpos(1 : archive.NP);
% 	  archive.pop = archive.pop(rndpos, :);
% 	end
%     end
    bsf_fit_var
%     conv(index,1) = nfes;
%     conv(index,2) = bsf_fit_var;
%     conv(index,3) = pop_size;
%     index = index+1;
%     end
%     conv(all(conv==0,2),:)=[];
%     bsf_error_val = bsf_fit_var;
%     if bsf_error_val < val_2_reach
%         bsf_error_val = 0;
%      end
% 
%     fprintf('%d th run, best-so-far error value = %1.8e\n', run_id , bsf_error_val)
%     
  end %% end 1 run
 
end %% end 1 function run
toc
Opt_value = bsf_solution;
Func(Opt_value);
Total_Turbine = size(wind_farm,1);
%wind_farm
fprintf('\n total_power %.6f',farm_power);
fprintf('\n Efficiency %.6f',Efficiency);
fprintf('\n Cost per kW %.8f',cost_per_kW);
%fprintf('\n x=%d   y=%d   r=%d   h=%d\n',wind_farm');
fprintf('\n Total Turbine %d',Total_Turbine);

%plot windfarm
    figure('DefaultAxesFontSize',9);
    x1 = wind_farm(:,1);
    y1 = wind_farm(:,2);
    plot_diameter = 2*wind_farm(:,3);
    plot_hubheight = wind_farm(:,4);
    plot(x1,y1,'m>');
    set(gca,'DefaultTextFontSize',9);
    hold on;
    text(x1,y1,num2str(plot_diameter),'VerticalAlignment','top');
    text(x1,y1,num2str(plot_hubheight),'VerticalAlignment','bottom');
    xlabel('x-coordinate','Fontsize',9);
    ylabel('y-coordinate','Fontsize',9);
    grid on;
