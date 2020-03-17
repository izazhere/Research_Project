'NS_ESPSO'
clear all
clc
n = 50;              % Size of the swarm " no of birds "
bird_step  = 10000;  % Maximum number of "birds steps" iterations 
maxfe = 200000;

dim = 30;            % Dimension of the problem state

ns_spsoplot = zeros(maxfe/5000 + 1,13);
cput = zeros(13,1);
succ = zeros(13*3, 1);
stat = zeros(13*4, 1);
runnum = 1;

for problem = 1 : 13
    funcid = problem;
  
  switch problem
        case 1
            lu = [-100 * ones(1, dim); 100 * ones(1, dim)];
            threshold = 1e-10;
        case 2
            lu = [-10 * ones(1, dim); 10 * ones(1, dim)];
            threshold = 1e-10;
        case 3
            lu = [-100 * ones(1, dim); 100 * ones(1, dim)];
            threshold = 1e-5;
        case 4
            lu = [-100 * ones(1, dim); 100 * ones(1, dim)];
            threshold = 1e-10;
        case 5
            lu = [-30* ones(1, dim); 30 * ones(1, dim)];
            threshold = 1e2;
        case 6
            lu = [-100 * ones(1, dim); 100 * ones(1, dim)];
             threshold = 1e-10;
        case 7
            lu = [-1.28 * ones(1, dim); 1.28 * ones(1, dim)];
            threshold = 0.1;
        case 8
            lu = [-500 * ones(1, dim); 500 * ones(1, dim)];
            threshold = 2e3;
        case 9
            lu = [-5.12 * ones(1, dim); 5.12 * ones(1, dim)];
            threshold = 1e2;
        case 10
            lu = [-32 * ones(1, dim); 32 * ones(1, dim)];
             threshold = 1e-10;
        case 11
            lu = [-600 * ones(1, dim); 600 * ones(1, dim)];
             threshold = 1e-10;

        case 12

            lu = [-50 * ones(1, dim); 50 * ones(1, dim)];
             threshold = 1e-10;
        case 13
            lu = [-50 * ones(1, dim); 50 * ones(1, dim)];
             threshold = 1e-10;
  end

    lu = lu';

        N_inteval=4;      %memebership interval
        c=zeros(2,N_inteval);
        
        %c(1,:)=[2,2.05,2.1,2.15,2.2,2.25,1.8,1.85];
        %c(2,:)=[2,1.95,1.9,1.85,1.8,1.75,2.2,2.15];
        
       c(1,:)=[2,2.1,2.2,1.8];
       c(2,:)=[2,1.9,1.8,2.2];
       w = 0.9;
       
    data = zeros(runnum,1);
    %several runs
for run = 1: runnum
    succornot = 0;
       
    % initialize parameters                                                               
        mv = 0.2 * (lu(:, 2) - lu(:, 1));
        VRmin = repmat(lu(:, 1), 1, n);
        VRmax = repmat(lu(:, 2), 1, n);
        Vmin = repmat(-mv, 1, n);
        Vmax = -Vmin;                         
        current_fitness = 0*ones(n,1);

        % initialize swarm and velocity 
                                 
        XRRmin = repmat(lu(:, 1), 1, n);
        XRRmax = repmat(lu(:, 2), 1, n);
        rand('seed', sum(100 * clock));
        current_position = [XRRmin + (XRRmax - XRRmin) .* rand(dim, n)];
        velocity = [XRRmin + (XRRmax - XRRmin) .* rand(dim, n)];
        local_best_position  = current_position ;

        % evaluate initial swarm

        current_fitness = test_func(current_position', problem);

        local_best_fitness  = current_fitness ;
        [global_best_fitness,g] = min(local_best_fitness) ;

        for i=1:n
            globl_best_position(:,i) = local_best_position(:,g) ;
        end
        
        % calculate evolutionary factor (E_f) and inertia (w)
        
        D_matix=squareform(pdist(current_position'));
        d=(sum(D_matix))'/n; 
        E_f=(d(g)-min(d))/(max(d)-min(d));
        %w=0.5*E_f+0.4;                  % pso momentum or inertia 

        % calculate states (ksi)on the bases of (E_f)
        
        for i=0:N_inteval                                
            if((i/N_inteval)<=E_f<((i+1)/N_inteval))
                ksi=i; 
            end
        end
                                               
        R1 = rand(dim, n);
        R2 = rand(dim, n);
    
        %velocity update equation

        velocity = w *velocity + c(1,ksi)*(R1.*(local_best_position-current_position)) + c(2,ksi)*(R2.*(globl_best_position-current_position));

        % position update equation 
                                                           
        current_position = current_position + velocity ;
                                          
        FES = n;
        gen = 0;
        ploti = 1;                                          
                
        %% Main Loop
        iter = 0 ;        % Iterations’counter
        tic;
        while  ( FES < maxfe )
        iter = iter + 1;

        if(mod(FES, 10000) == 0)
            w = 0.9 - (0.9 - 0.5)*(FES/maxfe);
            %w = 0.5*E_f+0.4*(FES/maxfe);
        end

        R1 = rand(dim, n);
        R2 = rand(dim, n);
        rand('seed', sum(100 * clock));

        current_fitness = test_func(current_position', problem);
        FES = FES + n;

        for i = 1 : n
                if current_fitness(i) < local_best_fitness(i)
                   local_best_fitness(i)  = current_fitness(i);  
                   local_best_position(:,i) = current_position(:,i)   ;
                end   
         end

         [current_global_best_fitness,g] = min(local_best_fitness);


        if current_global_best_fitness < global_best_fitness
           global_best_fitness = current_global_best_fitness;

            for i=1:n
                globl_best_position(:,i) = local_best_position(:,g);
            end

        end

        fprintf('Best fitness: %e\n', global_best_fitness); 

         velocity = w *velocity + c(1,ksi)*(R1.*(local_best_position-current_position)) + c(2,ksi)*(R2.*(globl_best_position-current_position));

         %velocity boundary
                for i = 1:n
                    velocity(:, i) = max(velocity(:,i), Vmin(:,i));
                    velocity(:, i) = min(velocity(:,i), Vmax(:,i));
                end

         current_position = current_position + velocity; 

         %position boundary
                for i = 1:n
                    current_position(:,i) = max(current_position(:,i), lu(:,1));
                    current_position(:,i) = min(current_position(:,i), lu(:,2));
                end

        %record plots
                if( ceil(FES/5000) == ploti )
                    ns_spsoplot(ploti, problem) = ns_spsoplot(ploti, funcid) + global_best_fitness;
                    ploti = ploti + 1;
                end;
         %threshold 
                if( global_best_fitness  <= threshold && ~succornot)
                    succornot = 1;
                    succ((funcid - 1)*3 + 1) = succ((funcid - 1)*3 + 1) + 1;
                    succ((funcid - 1)*3 + 2) = succ((funcid - 1)*3 + 2) + FES;
                end;

        x=current_position(1,:);
        y=current_position(2,:);

        end 
                     
       [Jbest_min,I] = min(current_fitness); % minimum fitness
                   current_position(:,I); % best solution

        fprintf('Run No.%d Done!\n', run); 
        disp(['CPU time: ',num2str(toc)]);
        cput(funcid) = cput(funcid) + toc;
        data(run) = global_best_fitness;
    
end

stat((funcid - 1)*4 + 1) = mean(data);
stat((funcid - 1)*4 + 2) = min(data);
stat((funcid - 1)*4 + 3) = std(data);

end 

save NS_ESPSO ns_spsoplot


% for i = 1:13
% figure;
% plot([0:5000:FES],[log10(ns_spsoplot(:,i)')],'gd-');hold on;
% end

xlswrite('D:\NS_ESPSO_Results\ns_spsoplot.xls',ns_spsoplot./runnum,'NSJPSO');
xlswrite('D:\NS_ESPSO_Results\time.xls',cput./runnum,'NSJPSO');
xlswrite('D:\NS_ESPSO_Results\succ.xls',succ./runnum,'NSJPSO');
xlswrite('D:\NS_ESPSO_Results\stat.xls',stat,'NSJPSO');        
         