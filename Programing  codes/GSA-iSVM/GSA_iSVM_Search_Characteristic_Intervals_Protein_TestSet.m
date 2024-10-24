%% a litte clean work
%close all;
clear;
clc;
format compact;
%% ��������
load nirs_rice_powder_data.mat;%153����������(2014������)
nirs_data=nirs_rice_powder_data;
content_data=xlsread('rice_protein_data.xls');  %153��label����(2014������)
%% ����Ԥ����
%SGƽ����Savitzky Golay filter��ƽ���˲�
% nirs_sg_data = savgol(nirs_data,7,3,0);
% nirs_data=nirs_sg_data;
%MSC��Ԫɢ��У����Multiplicative scatter correction��
  nirs_msc_data=msc(nirs_data,nirs_data);
  nirs_data=nirs_msc_data;
% %SNV��׼����任(Standard normal variate)
  nirs_snv_data=snv(nirs_data);
  nirs_data=nirs_snv_data;
% %���������ο�diff(x,n)����,�����������һ��
%    nirs_ds_data=diff(nirs_data',1);
%    nirs_data=nirs_ds_data';
%   nirs_ds_data=diff(nirs_data',1);
%   nirs_data=nirs_ds_data';
data=nirs_data;
%% �쳣�����޳�
outlier_num=[1 27 28 29 30 32 80]
data(outlier_num,:)=[];%ȥ���쳣������������
content_data(outlier_num,:)=[];%ȥ���쳣��������ֵ
%% ���ض������Լ����
%RS��:
% % % % Independent_validation_set_numbers=rs(data,23);
load Independent_validation_set_numbers_rice_powder.mat;
Independent_test_data=data(Independent_validation_set_numbers,:);   %��������������m��Ӧ����
independent_test_label=content_data(Independent_validation_set_numbers,:);
data(Independent_validation_set_numbers,:)=[];%ȥ���쳣������������
content_data(Independent_validation_set_numbers,:)=[];%ȥ���쳣��������ֵ
%% ����У������Ԥ�⼯
[m,dminmax]=ks(data,90); %mΪ��ѡ��������������
train_data=data(m,:);   %��������������m��Ӧ����
train_label=content_data(m,:);
data(m,:)=[]; %ɾ�������е�����m��Ӧ����
content_data(m,:)=[];
test_data=data;
test_label=content_data;
%% SVM����
TrainL = train_label;
Train = train_data;
TestL = [test_label;independent_test_label];
Test =[test_data;Independent_test_data];
%% ��һ��
[Train,Test] = scaleForSVM(Train,Test,-1,1);
[TrainL,TestL,ps] = scaleForSVM(TrainL,TestL,-1,1);
%% ��������
intervals_num=58;
intervals_index=intervals_divide(nirs_data,intervals_num);
%% GSA-iSVM������������
gsa_option.inittempFactor=100;
gsa_option.lowertempFactor=0.9;
gsa_option.maxgen =100;%����Ŵ�����
gsa_option.sizepop = intervals_num+60;%��Ⱥ����Ϊ�����������+SVM��������λ��
gsa_option.v = 10;
gsa_option.ggap = 0.9;
gsa_option.cbound = [0,100];
gsa_option.gbound = [0,100];
gsa_option.pbound = [0.001,1];    
RESULT=[];
count=1;
searchNum=2;%��������
while count<=searchNum
    [bestCVmse,bestc,bestg,bestp,BestWaveIndexs,BestIntervasString] = gsaSVMcgpANDintervlMultibitsForRegress(TrainL,Train,gsa_option,count,intervals_index); 
    RESULT=[RESULT; [intervals_num,length(BestWaveIndexs),length(BestIntervasString),bestCVmse,bestc,bestg,bestp]]; %��������ͳ����Ϣ
    for i=1:length(BestIntervasString)
        RESULTBestIntervasString(count,i)=BestIntervasString(i); 
    end
    count=count+1;    
end
RESULTBestIntervasString(RESULTBestIntervasString==0)=NaN; %������������

