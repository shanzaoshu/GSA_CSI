%% a litte clean work
%close all;
clear;
clc;
format compact;
%% 导入数据
load nirs_rice_powder_data.mat;%153个光谱数据(2014年数据)
nirs_data=nirs_rice_powder_data;
content_data=xlsread('rice_protein_data.xls');  %153个label数据(2014年数据)
%% 光谱预处理
%SG平滑：Savitzky Golay filter，平滑滤波
% nirs_sg_data = savgol(nirs_data,7,3,0);
% nirs_data=nirs_sg_data;
%MSC多元散射校正（Multiplicative scatter correction）
  nirs_msc_data=msc(nirs_data,nirs_data);
  nirs_data=nirs_msc_data;
% %SNV标准正则变换(Standard normal variate)
  nirs_snv_data=snv(nirs_data);
  nirs_data=nirs_snv_data;
% %导数处理，参考diff(x,n)函数,导数计算后少一列
%    nirs_ds_data=diff(nirs_data',1);
%    nirs_data=nirs_ds_data';
%   nirs_ds_data=diff(nirs_data',1);
%   nirs_data=nirs_ds_data';
data=nirs_data;
%% 异常样本剔除
outlier_num=[1 27 28 29 30 32 80]
data(outlier_num,:)=[];%去掉异常样本光谱数据
content_data(outlier_num,:)=[];%去掉异常样本含量值
%% 加载独立测试集编号
%RS法:
% % % % Independent_validation_set_numbers=rs(data,23);
load Independent_validation_set_numbers_rice_powder.mat;
Independent_test_data=data(Independent_validation_set_numbers,:);   %保留矩阵中数列m对应的行
independent_test_label=content_data(Independent_validation_set_numbers,:);
data(Independent_validation_set_numbers,:)=[];%去掉异常样本光谱数据
content_data(Independent_validation_set_numbers,:)=[];%去掉异常样本含量值
%% 构建校正集和预测集
[m,dminmax]=ks(data,90); %m为所选光谱行向量序列
train_data=data(m,:);   %保留矩阵中数列m对应的行
train_label=content_data(m,:);
data(m,:)=[]; %删除矩阵中的数列m对应的行
content_data(m,:)=[];
test_data=data;
test_label=content_data;
%% SVM数据
TrainL = train_label;
Train = train_data;
TestL = [test_label;independent_test_label];
Test =[test_data;Independent_test_data];
%% 归一化
[Train,Test] = scaleForSVM(Train,Test,-1,1);
[TrainL,TestL,ps] = scaleForSVM(TrainL,TestL,-1,1);
%% 划分区间
intervals_num=58;
intervals_index=intervals_divide(nirs_data,intervals_num);
%% GSA-iSVM搜索特征区间
gsa_option.inittempFactor=100;
gsa_option.lowertempFactor=0.9;
gsa_option.maxgen =100;%最大遗传代数
gsa_option.sizepop = intervals_num+60;%种群个数为划分区间个数+SVM参数编码位数
gsa_option.v = 10;
gsa_option.ggap = 0.9;
gsa_option.cbound = [0,100];
gsa_option.gbound = [0,100];
gsa_option.pbound = [0.001,1];    
RESULT=[];
count=1;
searchNum=2;%搜索次数
while count<=searchNum
    [bestCVmse,bestc,bestg,bestp,BestWaveIndexs,BestIntervasString] = gsaSVMcgpANDintervlMultibitsForRegress(TrainL,Train,gsa_option,count,intervals_index); 
    RESULT=[RESULT; [intervals_num,length(BestWaveIndexs),length(BestIntervasString),bestCVmse,bestc,bestg,bestp]]; %特征谱区统计信息
    for i=1:length(BestIntervasString)
        RESULTBestIntervasString(count,i)=BestIntervasString(i); 
    end
    count=count+1;    
end
RESULTBestIntervasString(RESULTBestIntervasString==0)=NaN; %特征谱区序列

