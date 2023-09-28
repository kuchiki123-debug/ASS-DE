function All_Run()
clc
clear all
close all
format long
format compact
addpath('Public');
addpath('DEs');
addpath('CEC2015');%Fun_Num=25;
TEV=Error(Dim(1));
runmax=25 ;
Fun_Num=1:15;

%算法内部参数
N=100;               %population size

%% ASS_DE
name = "ASS_DE";
for problem=1:15
    [RunResult,RunValue,RunTime,RunFES,RunOptimization,RunParameter]=ASS_DE(problem,N,runmax);
    sign=(RunResult<=TEV(problem));
    Ratio=sum(sign)/runmax;
    FES=sign.*RunFES;
    TestFitness=[TestFitness;RunResult];
    TestResult=[TestResult;min(RunResult) max(RunResult) median(RunResult) mean(RunResult) std(RunResult)];
    TestValue=[TestValue;mean(RunValue)];
    TestTime=[TestTime;mean(RunTime)];
    TestRatio=[TestRatio;Ratio];
    TestFES=[TestFES;mean(FES)];
    TestOptimization=[TestOptimization;RunOptimization];
    TestParameter=[TestParameter;RunParameter];
end
Test=sprintf('TestFitness/%s.mat',name);
save(Test,'TestFitness');
Test=sprintf('TestResult/%s.mat',name);
save(Test,'TestResult');
Test=sprintf('TestValue_FES/%s.mat',name);
save(Test,'TestValue');
Test=sprintf('TestTime/%s.mat',name);
save(Test,'TestTime');
Test=sprintf('TestRatio/%s.mat',name);
save(Test,'TestRatio');
Test=sprintf('TestFES/%s.mat',name);
save(Test,'TestFES');
Test=sprintf('TestOptimization/%s.mat',name);
save(Test,'TestOptimization');
Test=sprintf('TestParameter_FES/%s.mat',name);
save(Test,'TestParameter');
end

