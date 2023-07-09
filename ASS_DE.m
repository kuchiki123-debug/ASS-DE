%% The original Different Evolution
%problem: the serial number of testing function recorded in "Public\benchmark_func.m"
%N: the population size
%runmax: the number of the algorithm runs
%RunResult: the  optimal value produced by each algorithm runs
%RunOptimization: the optimal value produced by reach algorithm runs
%RunValue: the fitness of optimal value produced by each 10000 FES
%RunParameter:the optimal value produced by each 10000 FES
%RunTime: the time spent by each algorithm runs
%RunFES: the FES required to satisfy the conditions
function [RunResult,RunValue,RunTime,RunFES,RunOptimization,RunParameter]=ASS_DE(img,imn,problem,N,runmax)
    'ASS_DE'
    D=Dim(problem);
    lu=Boundary(D);
    tempTEV=Error(D);
    TEV = tempTEV(problem);
    G=50;
    FESMAX=G*N;
    RunOptimization=zeros(runmax,D);
    for run=1:runmax
        TimeFlag=0;
        TempFES=FESMAX;
        t1=clock;
        %--------------------CR--------------------
        CR=zeros(N,1);%交叉概率
        good_uCR=0.5*ones(1,4);
        goodCR=zeros(N,4);
        %--------------------F---------------------
        F=zeros(N,1);%变异步长（变异尺度）
        good_uF=0.5*ones(1,4);
        goodF=zeros(N,4);
        %------------------------------------------
        gama=0.95;
        alpha=0.9;
        c=0.1;
        x=Initpop(N,D,lu);%种群初始化，参考CEP
        x_old=[];
        v=zeros(N,D);
        fitness=benchmark_func(imn,img,x);%计算每一个个体的函数值，参考CEP
        [x_best,x_best_number]=min(fitness);
%         min_fitness=x_best;
        stagnation=0;
        FES=N;%当前的函数评价次数，即函数已计算的次数
        k=1;
        p=0.05;
        %初始化每个个体的参数Q表
        Q=zeros(4,4,N);
        state=ceil(rand(N,1)*3);
        state1=0;
        state2=0;
        state3=0;
        state4=0;
        P=0;
        Flag=0;
        epsilon=0.5;
        while FES<=FESMAX
            Fitness = fitness;
            min_fitness=x_best;
            epsilon=epsilon-(1/(FESMAX/100*2));
            for i=1:N
               %% calculate the action selection probability
                pos=0;
                pso=0;
                if Flag==1
                    a=ceil(rand*4);
                elseif Flag==2
                    a=ceil(rand*3);
                elseif max(Q(state(i),:,i))==0
                    a=ceil(rand*4);
                elseif epsilon < rand
                    [~,a] = max(Q(state(i),:,i));
                else
                    a=ceil(rand*4);
                end
                
                % Generate CR according to a normal distribution with mean CRm, and std 0.1
                % Generate F according to a cauchy distribution with location parameter Fm, and scale parameter 0.1
                [F(i,a), CR(i,a)] = randFCR(good_uCR(1,a), 0.1, good_uF(1,a), 0.1);
                
               %%  mution %%%%
                x_all=[x; x_old];
                indexSet=1:size(x_all,1);%生成一个1,2,3，。。。N的序列
                indexSet(i)=[];%去掉向量indexSet第i个元素，即把i删除
                temp=floor(rand*(N-1))+1;%在1到N-1之间随机取一个整数
                r(1)=indexSet(temp);%生成差分变异的第一个个体角标
                indexSet(temp)=[];%去掉向量indexSet第temp个元素
                temp=floor(rand*(N-2))+1;
                r(2)=indexSet(temp);
                indexSet(temp)=[];
                temp=floor(rand*(N-3))+1;
                r(3)=indexSet(temp);
                indexSet(temp)=[];
                temp=floor(rand*(N-4))+1;
                r(4)=indexSet(temp);
                indexSet(temp)=[];
                temp=floor(rand*(N-5))+1;
                r(5)=indexSet(temp);
                indexSet(temp)=[];
                temp=floor(rand*(size(x_all,1)-6))+1;
                r(6)=indexSet(temp);
                % ---------------Find the p-best solutions-----------------
                pNP = max(round(p * N), 2); % choose at least two best solutions
                randindex=ceil(rand*pNP);
                [valBest, indBest] = sort(fitness, 'ascend');
                pbest = x(indBest(randindex), :); % randomly choose one of the top 100p% solutions
                if a==1
                    v(i,:)=x(i,:)+F(i,a)*(pbest-x(i,:))+F(i,a)*(x(r(1),:)-x_all(r(6),:));%DE/current-to-pbest/1
                elseif a==2
                    v(i,:)=x(i,:)+F(i,a)*(x(indBest(1),:)-x(i,:))+F(i,a)*(x(r(1),:)-x(r(2),:));%DE/current-to-best/1
                elseif a==3
                    v(i,:)=x(i,:)+F(i,a)*(x(r(1),:)-x(i,:))+F(i,a)*(x(r(2),:)-x_all(r(3),:));%DE/current-to-rand/1
                elseif a==4
                    v(i,:)=x(indBest(1),:)+F(i,a)*(x(r(2),:)-x(r(3),:));%best/1
                end

                
                %%  crossover %%%%
                jrand=floor(rand*D)+1;%在1-D之间随机取一维
                for j=1:D
                    if  v(i,j)>lu(2,j)
                        v(i,j)=max(lu(1,j),2*lu(2,j)-v(i,j));%超出上界处理，参考CEP
                    end
                    if  v(i,j)<lu(1,j)
                        v(i,j)=min(lu(2,j),2*lu(1,j)-v(i,j));%超出下界处理，参考CEP
                    end
                    if (rand<CR(i))||(j==jrand)
                        newx(i,j)=v(i,j);%二项式交叉
                    else
                        newx(i,j)=x(i,j);
                    end
                end
                newx(i,5) = round(newx(i,5));
                newx(i,6) = round(newx(i,6));
                children_fitness(i)=benchmark_func(imn,img,newx(i,:));%计算生成的试验向量newx(i,:)的目标函数值
 
                %% selection %%%%
                if children_fitness(i)<fitness(i)%贪婪选择
                    x(i,:)=newx(i,:);%个体更新  
                    promote_rate=(fitness(i)-children_fitness(i))/fitness(i);
                    P=[P,promote_rate];
                    fitness(i)=children_fitness(i);%函数值更新
                    if promote_rate<0.02
                        new_state=1;
                        state1=state1+1;
                        reward=promote_rate*100;
                    elseif promote_rate<0.05
                        new_state=2;
                        state2=state2+1;
                        reward=promote_rate*100;
                    elseif promote_rate<0.1
                        new_state=3;
                        state3=state3+1;
                        reward=promote_rate*100;
                    else
                        new_state=4;
                        state4=state4+1;
                        reward=promote_rate*100;
                    end
                    goodF(i,a) = F(i,a);
                    goodCR(i,a) = CR(i,a);
                    if children_fitness(i)<x_best
                        x(x_best_number,:)=newx(i,:);
                        x_best=children_fitness(i);
                    end
                else
                    reward=0;
                    new_state=ceil(rand*4);
                end          
        
                %% Update Q value
                Q(state(i,:),a,i)=Q(state(i,:),a,i)+alpha*(reward+gama*max(Q(new_state,:,i))-Q(state(i,:),a,i));
                
                state(i,:)=new_state;

                FES=FES+1;
                if mod(FES,50) == 0
                    [kk,ll]=min(fitness);
                    RunValue(run,k)=kk;
                    Para(k,:)=x(ll,:);
                    k=k+1;
                    fprintf('Algorithm:%s problemIndex:%d Run:%d FES:%d Best:%g\n','ASS_DE',problem,run,FES,kk);
                end
                if TimeFlag==0
                    if min(fitness)<=TEV
                        TempFES=FES;
                        TimeFlag=1;
                    end
                end
            end
            %----------------------------------F--------------------------------------
            for K=1:4
                if sum(goodF(:,K)) > 0
                    goodF_pos = find(goodF(:,K) > 0);
                    good_uF(1,K)=(1-c)*good_uF(1,K)+c*(sum((goodF(goodF_pos,K)).^2)/sum(goodF(goodF_pos,K)));
                end
            end
            %----------------------------------CR-------------------------------------
            for K=1:4
                if sum(goodCR(:,K)) > 0
                    goodCR_pos = goodCR(:,K) > 0;
                    good_uCR(1,K)=(1-c)*good_uCR(1,K)+c*(mean(goodCR(goodCR_pos,K)));
                end
            end
           
            for o=1:N
                if Fitness(o) > fitness(o)
                    pos=[pos,o];
                else
                    pso=[pso,o];
                end
            end

            if ~isempty(pso)
                pso(1)=[];
            end
            
            x_old=x(pso,:);
            
            goodF=zeros(N,4);
            goodCR=zeros(N,4);

        end
        [kk,ll]=min(fitness);
        gbest=x(ll,:);
        t2=clock;
        RunTime(run)=etime(t2,t1);
        RunResult(run)=kk;
        RunFES(run)=TempFES;
        RunOptimization(run,1:D)=gbest;
        RunParameter{run}=Para;
    end
end

function [F,CR] = randFCR( CRmgood_uCR, CRsigma, good_uF,  Fsigma)

% this function generate CR according to a normal distribution with mean "CRm" and sigma "CRsigma"
%           If CR > 1, set CR = 1. If CR < 0, set CR = 0.
% this function generate F  according to a cauchy distribution with location parameter "Fm" and scale parameter "Fsigma"
%           If F > 1, set F = 1. If F <= 0, regenrate F.
%
% Version: 1.1   Date: 11/20/2007
% Written by Jingqiao Zhang (jingqiao@gmail.com)

    %% generate CR
    CR = CRmgood_uCR + CRsigma * randn();
    CR = min(1, max(0, CR));                % truncated to [0 1]

    %% generate F
    F = randCauchy(good_uF, Fsigma);
    F = min(1, F);                          % truncation

    % we don't want F = 0. So, if F<=0, we regenerate F (instead of trucating it to 0)

    while F<=0
        F = randCauchy(good_uF, Fsigma);
        F = min(1, F);                      % truncation
    end
end

% Cauchy distribution: cauchypdf = @(x, mu, delta) 1/pi*delta./((x-mu).^2+delta^2)
function result = randCauchy(mu, delta)

% http://en.wikipedia.org/wiki/Cauchy_distribution
result = mu + delta * tan(pi * (rand - 0.5));
end