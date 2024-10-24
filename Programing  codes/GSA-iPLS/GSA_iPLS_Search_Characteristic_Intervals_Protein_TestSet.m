%% a litte clean work
%close all;
clear;
clc;
format compact;
%% ��������
load nirs_rice_powder_data.mat;%150����������
nirs_data=nirs_rice_powder_data;
content_data=xlsread('rice_protein_data.xls');  %150��label����
%% ����Ԥ����
%MSC��Ԫɢ��У����Multiplicative scatter correction��
  nirs_msc_data=msc(nirs_data,nirs_data);
  nirs_data=nirs_msc_data;
% %SNV��׼����任(Standard normal variate)
  nirs_snv_data=snv(nirs_data);
  nirs_data=nirs_snv_data;
data=nirs_data;
%% �쳣�����޳�
outlier_num=[1 27 28 29 30 32 80]
data(outlier_num,:)=[];%ȥ���쳣������������
content_data(outlier_num,:)=[];%ȥ���쳣��������ֵ
%RS��:
% Independent_validation_set_numbers=rs(data,23);
load Independent_validation_set_numbers_rice_powder.mat;
Independent_test_data=data(Independent_validation_set_numbers,:);   %��������������m��Ӧ����
independent_test_label=content_data(Independent_validation_set_numbers,:);
data(Independent_validation_set_numbers,:)=[];%ȥ���쳣������������
content_data(Independent_validation_set_numbers,:)=[];%ȥ���쳣��������ֵ
%% ����У������Ԥ�⼯
%KS��
[m,dminmax]=ks(data,90); %mΪ��ѡ��������������
train_data=data(m,:);   %��������������m��Ӧ����
train_label=content_data(m,:);
data(m,:)=[]; %ɾ�������е�����m��Ӧ����
content_data(m,:)=[];
test_data=data;
test_label=content_data;
%% ��������
intervals_num=59;
intervals_index=intervals_divide(nirs_data,intervals_num);
%% GSA-iPLS������������
gsa_option.inittempFactor=100;
gsa_option.lowertempFactor=0.90;%����ϵ��
gsa_option.maxgen =100;%����Ŵ�����
gsa_option.sizepop =intervals_num;%��Ⱥ��ģ
gsa_option.m=ceil(intervals_num/10); %�������Ŷ�λ��,����10��ȡ������ȡ֤��floor()
gsa_option.v = 10;%K�۽�����֤��K
gsa_option.ggap = 0.9;%ѡ������Ĵ���
RESULT=[];
count=1;
TotalSearchCount=2;%�����ܴ���
while count<=TotalSearchCount
    [BestWaveIndexs,BestIntervasString,RMSECV]=gsaPLSintervalsForRegress(train_data,train_label,gsa_option,intervals_index,count);
    RESULT=[RESULT; [intervals_num,length(BestWaveIndexs),length(BestIntervasString),RMSECV]];
    for i=1:length(BestIntervasString)
        RESULTBestIntervasString(count,i)=BestIntervasString(i);  %��������ͳ����Ϣ
    end
    count=count+1;
end
RESULTBestIntervasString(RESULTBestIntervasString==0)=NaN; %������������