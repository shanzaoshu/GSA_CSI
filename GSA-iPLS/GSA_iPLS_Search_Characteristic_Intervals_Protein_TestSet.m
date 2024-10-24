%% a litte clean work
%close all;
clear;
clc;
format compact;
%% 导入数据
load nirs_rice_powder_data.mat;%150个光谱数据
nirs_data=nirs_rice_powder_data;
content_data=xlsread('rice_protein_data.xls');  %150个label数据
%% 光谱预处理
%MSC多元散射校正（Multiplicative scatter correction）
  nirs_msc_data=msc(nirs_data,nirs_data);
  nirs_data=nirs_msc_data;
% %SNV标准正则变换(Standard normal variate)
  nirs_snv_data=snv(nirs_data);
  nirs_data=nirs_snv_data;
data=nirs_data;
%% 异常样本剔除
outlier_num=[1 27 28 29 30 32 80]
data(outlier_num,:)=[];%去掉异常样本光谱数据
content_data(outlier_num,:)=[];%去掉异常样本含量值
%RS法:
% Independent_validation_set_numbers=rs(data,23);
load Independent_validation_set_numbers_rice_powder.mat;
Independent_test_data=data(Independent_validation_set_numbers,:);   %保留矩阵中数列m对应的行
independent_test_label=content_data(Independent_validation_set_numbers,:);
data(Independent_validation_set_numbers,:)=[];%去掉异常样本光谱数据
content_data(Independent_validation_set_numbers,:)=[];%去掉异常样本含量值
%% 构建校正集和预测集
%KS法
[m,dminmax]=ks(data,90); %m为所选光谱行向量序列
train_data=data(m,:);   %保留矩阵中数列m对应的行
train_label=content_data(m,:);
data(m,:)=[]; %删除矩阵中的数列m对应的行
content_data(m,:)=[];
test_data=data;
test_label=content_data;
%% 划分区间
intervals_num=59;
intervals_index=intervals_divide(nirs_data,intervals_num);
%% GSA-iPLS搜索特征区间
gsa_option.inittempFactor=100;
gsa_option.lowertempFactor=0.90;%降温系数
gsa_option.maxgen =100;%最大遗传代数
gsa_option.sizepop =intervals_num;%种群规模
gsa_option.m=ceil(intervals_num/10); %邻域解的扰动位数,除以10上取整，下取证用floor()
gsa_option.v = 10;%K折交叉验证的K
gsa_option.ggap = 0.9;%选择操作的代沟
RESULT=[];
count=1;
TotalSearchCount=2;%搜索总次数
while count<=TotalSearchCount
    [BestWaveIndexs,BestIntervasString,RMSECV]=gsaPLSintervalsForRegress(train_data,train_label,gsa_option,intervals_index,count);
    RESULT=[RESULT; [intervals_num,length(BestWaveIndexs),length(BestIntervasString),RMSECV]];
    for i=1:length(BestIntervasString)
        RESULTBestIntervasString(count,i)=BestIntervasString(i);  %特征谱区统计信息
    end
    count=count+1;
end
RESULTBestIntervasString(RESULTBestIntervasString==0)=NaN; %特征谱区序列