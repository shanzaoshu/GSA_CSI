function [BestMSE,Bestc,Bestg,Bestp,BestWaveIndexs,BestIntervalsString] = gsaSVMcgpANDintervlMultibitsForRegress(train_label,train_data,gsa_option,k,intervals_index)
%%SVM参数和特征波长同步优化
% gsaSVMcpgForClass
%
% by faruto
%Email:patrick.lee@foxmail.com QQ:516667408 http://blog.sina.com.cn/faruto
%last modified 2011.06.08
% 若转载请注明：
% faruto and liyang , LIBSVM-farutoUltimateVersion 
% a toolbox with implements for support vector machines based on libsvm, 2011. 
% Software available at http://www.matlabsky.com23 
% 
% Chih-Chung Chang and Chih-Jen Lin, LIBSVM : a library for
% support vector machines, 2001. Software available at
% http://www.csie.ntu.edu.tw/~cjlin/libsvm
%% 参数初始化
if nargin == 2
    gsa_option = struct('inittempFactor',100,'lowertempFactor',0.8,'maxgen',200,'sizepop',20,'ggap',0.9,...
        'cbound',[0,100],'gbound',[0,1000],'pbound',[0.01,1],'v',5);
end
% inittemp:初始温度，默认为100
%lowertemp:降温系数，默认为0.8
% maxgen:最大的进化代数,默认为200,一般取值范围为[100,500]
% sizepop:种群最大数量,默认为20,一般取值范围为[20,100]
% cbound = [cmin,cmax],参数c的变化范围,默认为(0,100]
% gbound = [gmin,gmax],参数g的变化范围,默认为[0,1000]
% pbound = [pmin,pmax],参数p的变化范围,默认为[0,1]
% v:SVM Cross Validation参数,默认为5

%%
MAXGEN = gsa_option.maxgen;
NIND = gsa_option.sizepop;
NVAR = 3;  %特征区间独立编码，其他编码为60位染色体
PRECI = 20;
IntervalCodeLength=size(intervals_index,1);
GGAP = gsa_option.ggap;
trace = zeros(MAXGEN,2);

FieldID = ...
[rep([PRECI],[1,NVAR]); ...
[gsa_option.cbound(1),gsa_option.gbound(1),gsa_option.pbound(1);gsa_option.cbound(2),gsa_option.gbound(2),gsa_option.pbound(2);];...
[0,0,0;0,0,0;0,1,1;1,1,1]];

%Create an initial population
Chrom = crtbp(NIND,IntervalCodeLength+NVAR*PRECI);%码长为特征区间长度+60,染色体前面为特征区间，后面为60位支持向量机参数
gen = 1;
v = gsa_option.v;
BestMSE = inf;
Bestc = 0;
Bestg = 0;
Bestp = 0;
%% Binary string to real vector
cg = bs2rv(Chrom(:,IntervalCodeLength+1:end),FieldID);

for nind = 1:NIND
    cmd = ['-v ',num2str(v),' -c ',num2str(cg(nind,1)),' -g ',num2str(cg(nind,2)),' -p ',num2str(cg(nind,3)),' -s 3'];
    %求bestwaveindexs_temp
    bestwaveindexs_temp=[];
    intervasString=Chrom(nind,1:IntervalCodeLength);
    for i=1:IntervalCodeLength
        if intervasString(i)==1
            for j=intervals_index(i,1):intervals_index(i,2)
               bestwaveindexs_temp=[bestwaveindexs_temp j];
            end
        end
    end
    train_data_best=train_data(:,bestwaveindexs_temp);
    %ObjV为列向量，适应度函数值
    ObjV(nind,1) = svmtrain(train_label,train_data_best,cmd);
end
%输出最优解及其序号，Y为最优解，I为染色体在种群中的序号
[BestMSE,I] = min(ObjV);
[BadMSE,J]=max(ObjV);
Bestc = cg(I,1);
Bestg = cg(I,2);
Bestp = cg(I,3);
BestWaveIndexs=[];
BestIntervalsString=[];
intervasStringbest=Chrom(I,1:IntervalCodeLength);
for i=1:IntervalCodeLength
    if intervasStringbest(i)==1
        BestIntervalsString=[BestIntervalsString i];
        for j=intervals_index(i,1):intervals_index(i,2)
           BestWaveIndexs=[BestWaveIndexs j];
        end
    end
end

%模拟退火初温确定
INITTEMP=gsa_option.inittempFactor*(BadMSE-BestMSE);
%退温系数
LOWERTEMP=gsa_option.lowertempFactor;
% disp(INITTEMP);
% disp(LOWERTEMP);

steps=MAXGEN;
hwait=waitbar(0,'请等待>>>>>>>>');%进度条
step=steps/100;
while 1
	%% 进度条
    if steps-gen<=1
        str1=['第', num2str(k),'次','即将完成'];
        waitbar(gen/steps,hwait,str1);
        pause(0.05);
    else
        PerStr=fix(gen/step);
        str2=['第', num2str(k),'次','正在进行中',num2str(PerStr),'%'];
        waitbar(gen/steps,hwait,str2);
        pause(0.05);
    end

    %个体适应度值计算，结合温度参数的适应度函数设计
    for chromNum = 1:size(ObjV,1)
        FitnV(chromNum,1)=exp(-(ObjV(chromNum,1)-min(ObjV))/INITTEMP);    
    end
    %遗传操作
    SelCh = select('rws',Chrom,FitnV,GGAP);
    SelCh = recombin('xovsp',SelCh,0.7);
    SelCh = mut(SelCh);
    %模拟退火操作：特征区间多维变异
    SelCh = simulatedAnneal_intervals_multibits(SelCh,INITTEMP,FieldID,train_label,train_data,v,intervals_index);  
    
 %新子代目标函数MSE计算 
    cg = bs2rv(SelCh(:,IntervalCodeLength+1:end),FieldID);
    for nind = 1:size(SelCh,1)
        cmd = ['-v ',num2str(v),' -c ',num2str(cg(nind,1)),' -g ',num2str(cg(nind,2)),' -p ',num2str(cg(nind,3)),' -s 3'];
        %求bestWaveIndexs_temp
        bestWaveIndexs_temp=[];
        intervasString=SelCh(nind,1:IntervalCodeLength);
        for i=1:IntervalCodeLength
            if intervasString(i)==1
                for j=intervals_index(i,1):intervals_index(i,2)
                   bestWaveIndexs_temp=[bestWaveIndexs_temp j];
                end
            end
        end
        train_data_best=train_data(:,bestwaveindexs_temp);     
        
        ObjVSel(nind,1) = svmtrain(train_label,train_data_best,cmd);
    end
    %将子代插入上一代种群
    [Chrom,ObjV] = reins(Chrom,SelCh,1,1,ObjV,ObjVSel);   
    
    [NewBestCVaccuracy,I] = min(ObjV);
    cg_temp = bs2rv(Chrom(:,IntervalCodeLength+1:end),FieldID);
    temp_NewBestCVaccuracy = NewBestCVaccuracy;
    
    if NewBestCVaccuracy < BestMSE
       BestMSE = NewBestCVaccuracy;
       Bestc = cg_temp(I,1);
       Bestg = cg_temp(I,2);
       Bestp = cg_temp(I,3);
       BestWaveIndexs=[];
       BestIntervalsString=[];
       intervasStringbest=Chrom(I,1:IntervalCodeLength);
        for i=1:IntervalCodeLength
            if intervasStringbest(i)==1
                BestIntervalsString=[BestIntervalsString i];
                for j=intervals_index(i,1):intervals_index(i,2)
                   BestWaveIndexs=[BestWaveIndexs j];
                end
            end
        end
    end
  %尽量使svm的参数c最小  
    if abs( NewBestCVaccuracy-BestMSE ) <= 10^(-4) && ...
        cg_temp(I,1) < Bestc
       BestMSE = NewBestCVaccuracy;
       Bestc = cg_temp(I,1);
       Bestg = cg_temp(I,2);
       Bestp = cg_temp(I,3);
       BestWaveIndexs=[];
       BestIntervalsString=[];
       intervasStringbest=Chrom(I,1:IntervalCodeLength);
        for i=1:IntervalCodeLength
            if intervasStringbest(i)==1
                BestIntervalsString=[BestIntervalsString i];
                for j=intervals_index(i,1):intervals_index(i,2)
                   BestWaveIndexs=[BestWaveIndexs j];
                end
            end
        end
    end    
    %当前带最佳适应度和平均适应度
    trace(gen,1) = min(ObjV);
    trace(gen,2) = sum(ObjV)/length(ObjV);
    
%     if gen >= MAXGEN/2 && ...
%        ( temp_NewBestCVaccuracy-BestMSE ) <= 10^(-4)
%         break;
%     end
    if gen == MAXGEN
        break;
    end
    gen = gen + 1;
    INITTEMP=LOWERTEMP*INITTEMP;
end
close(hwait);%关闭进度条
