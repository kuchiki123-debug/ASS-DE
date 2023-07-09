function All_Run()
clc
clear all
close all
format long
format compact
addpath('Public');
addpath('denoising');
addpath('Test image');
TEV=Error(Dim(1));
runmax=5;
Fun_Num=2;
%Read image and add noise
imn = zeros(550,1000,3,3);
img = (im2double((imread('路径：****'))));%im2double其实就是double(I/255);像素值被标准化到0—1。（更换图片对应路径即可）
imn1 = imnoise(img,'salt & pepper',0.2);
imn2 = imnoise(img,'gaussian',0,0.2);
imn3 = imnoise(img,'speckle',0.2);
imn(:,:,:,1) = imn1;imn(:,:,:,2) = imn2;imn(:,:,:,3) = imn3;

%算法内部参数
N=10;%population size
Datime = date;
TestFitness=[];
TestResult=[];
TestValue={};
TestTime=[];
TestRatio=[];
TestFES=[];
TestOptimization={};
TestParameter={};
%% DE
name = "ASS_DE";
for problem=1:Fun_Num
    [RunResult,RunValue,RunTime,RunFES,RunOptimization,RunParameter]=ASS_DE(img,imn(:,:,:,problem),problem,N,runmax);
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


