function [BestWavelengthIndexs,BestIntervalsString,BestMSE] = gsaPLSintervalsForRegress(train_data,train_label,gsa_option,intervals_index,k)
%% 参数初始化
%train_label为某一单一成分值，以该值为基准进行优化
if nargin == 2
    gsa_option = struct('inittempFactor',100,'lowertempFactor',0.8,'maxgen',200,'sizepop',30,'ggap',0.9,'m',10,'v',5);
end
% inittemp:初始温度，默认为100
%lowertemp:降温系数，默认为0.8
% maxgen:最大的进化代数,默认为200,一般取值范围为[100,500]
% sizepop:种群最大数量,默认为20,一般取值范围为[20,100]
% v:SVM Cross Validation参数,默认为5

MAXGEN = gsa_option.maxgen;
NIND = gsa_option.sizepop;
IntervalCodeLength=size(intervals_index,1);
GGAP = gsa_option.ggap;
trace = zeros(MAXGEN,4);%依次保存最佳目标函数、平均目标函数、最大适应度函数和最小适应度函数值
v = gsa_option.v;
m =gsa_option.m;

%% Create an initial population
Chrom = crtbp(NIND,IntervalCodeLength);


gen = 1;
BestMSE = inf; %inf表示无穷大
BadMSE =inf;

%% 计算种群目标函数值
for nind = 1:NIND
    %ObjV为列向量，适应度函数值RMSECV,要计算每一个染色体的RMSE
    bestwaveindexs_temp=[];
    intervasString=Chrom(nind,:);
    for i=1:IntervalCodeLength
        if intervasString(i)==1
            for j=intervals_index(i,1):intervals_index(i,2)
               bestwaveindexs_temp=[bestwaveindexs_temp j];
            end
        end
    end
    ObjV(nind,1) = pls_crossvalind_wavelength(train_data,train_label,bestwaveindexs_temp,v,10);
end
%输出最优解及其序号，Y为最优解，I为染色体在种群中的序号
[BestMSE,I] = min(ObjV);
[BadMSE,J]=max(ObjV);
BestIntervalsString=[];
BestWavelengthIndexs=[];
intervasStringbest=Chrom(I,:);
for i=1:IntervalCodeLength
    if intervasStringbest(i)==1
        BestIntervalsString=[BestIntervalsString i];
        for j=intervals_index(i,1):intervals_index(i,2)
           BestWavelengthIndexs=[BestWavelengthIndexs j];
        end
    end
end
%模拟退火初温确定
INITTEMP=gsa_option.inittempFactor*(BadMSE-BestMSE);
%退温系数
LOWERTEMP=gsa_option.lowertempFactor;
%% 
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
  
    %结合温度参数的适应度函数设计
    for chromNum = 1:size(ObjV,1)
        FitnV(chromNum,1)=exp(-(ObjV(chromNum,1)-min(ObjV))/INITTEMP);    
    end
    %遗传操作
    SelCh = select('rws',Chrom,FitnV,GGAP);%赌轮选择
   % SelCh = recombin('xovsp',SelCh,0.7); %单点交叉
    SelCh = recombin('recdis',SelCh,0.7); %离散重组
    SelCh = mut(SelCh,0.01);%离散变异
    %模拟退火操作
    SelCh = pls_interval_simulatedAnneal(SelCh,INITTEMP,train_data,train_label,m,intervals_index,v);  
    
 %新子代目标函数MSE计算 
    for nind = 1:size(SelCh,1)
        bestwaveindexs_temp=[];
        intervasString=SelCh(nind,:);
        for i=1:IntervalCodeLength
            if intervasString(i)==1
                for j=intervals_index(i,1):intervals_index(i,2)
                   bestwaveindexs_temp=[bestwaveindexs_temp j];
                end
            end
        end
        ObjVSel(nind,1) = pls_crossvalind_wavelength(train_data,train_label,bestwaveindexs_temp,v,10);
       
    end
    %将子代插入上一代种群
    [Chrom,ObjV] = reins(Chrom,SelCh,1,1,ObjV,ObjVSel);   
    
    [NewBestCVaccuracy,I] = min(ObjV);
    temp_NewBestCVaccuracy = NewBestCVaccuracy;
    
    if NewBestCVaccuracy < BestMSE
           BestMSE = NewBestCVaccuracy;
           BestWavelengthIndexs=[];
           BestIntervalsString=[];
           intervasStringbest=Chrom(I,:);
            for i=1:IntervalCodeLength
                if intervasStringbest(i)==1
                    BestIntervalsString=[BestIntervalsString i];
                    for j=intervals_index(i,1):intervals_index(i,2)
                       BestWavelengthIndexs=[BestWavelengthIndexs j];
                    end
                end
            end
    end
      %当前带最佳适应度和平均适应度
    trace(gen,1) = min(ObjV);
    trace(gen,2) = sum(ObjV)/length(ObjV);
    trace(gen,3) = sum(FitnV)/length(FitnV);
    trace(gen,4) = min(FitnV);
    
    if gen == MAXGEN
        break;
    end
    gen = gen + 1;
    INITTEMP=LOWERTEMP*INITTEMP;
	
end

  close(hwait);%关闭进度条
%%
% figure;
% hold on;
% set(gcf,'color','white');
% 
% %设置绘图大小
% set(gcf,'Units','centimeters','Position',[10 10 17 13]);
% 
% trace = round(trace*10000)/10000;
% plot(trace(1:gen,1),'k*-');
% plot(trace(1:gen,2),'kd-');
% 
% tl=legend('最佳目标函数值','平均目标函数值');
% %设置tl的字体大小
% set(tl,'fontsize',9);
% xlabel('进化代数','FontSize',9);
% ylabel('目标函数值','FontSize',9);
% %设置坐标轴刻度值字体大小
% 
% figure;
% hold on;
% set(gcf,'color','white');
% 
% %设置绘图大小
% set(gcf,'Units','centimeters','Position',[10 10 17 13]);
% 
% trace = round(trace*10000)/10000;
% plot(trace(1:gen,3),'k*-');
% plot(trace(1:gen,4),'kd-');
% 
% tl=legend('平均适应度','最小适应度');
% %设置tl的字体大小
% set(tl,'fontsize',9);
% xlabel('进化代数','FontSize',9);
% ylabel('适应度函数值','FontSize',9);
% %设置坐标轴刻度值字体大小












