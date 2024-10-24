function [BestMSE,Bestc,Bestg,Bestp,BestWaveIndexs,BestIntervalsString] = gsaSVMcgpANDintervlMultibitsForRegress(train_label,train_data,gsa_option,k,intervals_index)
%%SVM��������������ͬ���Ż�
% gsaSVMcpgForClass
%
% by faruto
%Email:patrick.lee@foxmail.com QQ:516667408 http://blog.sina.com.cn/faruto
%last modified 2011.06.08
% ��ת����ע����
% faruto and liyang , LIBSVM-farutoUltimateVersion 
% a toolbox with implements for support vector machines based on libsvm, 2011. 
% Software available at http://www.matlabsky.com23 
% 
% Chih-Chung Chang and Chih-Jen Lin, LIBSVM : a library for
% support vector machines, 2001. Software available at
% http://www.csie.ntu.edu.tw/~cjlin/libsvm
%% ������ʼ��
if nargin == 2
    gsa_option = struct('inittempFactor',100,'lowertempFactor',0.8,'maxgen',200,'sizepop',20,'ggap',0.9,...
        'cbound',[0,100],'gbound',[0,1000],'pbound',[0.01,1],'v',5);
end
% inittemp:��ʼ�¶ȣ�Ĭ��Ϊ100
%lowertemp:����ϵ����Ĭ��Ϊ0.8
% maxgen:���Ľ�������,Ĭ��Ϊ200,һ��ȡֵ��ΧΪ[100,500]
% sizepop:��Ⱥ�������,Ĭ��Ϊ20,һ��ȡֵ��ΧΪ[20,100]
% cbound = [cmin,cmax],����c�ı仯��Χ,Ĭ��Ϊ(0,100]
% gbound = [gmin,gmax],����g�ı仯��Χ,Ĭ��Ϊ[0,1000]
% pbound = [pmin,pmax],����p�ı仯��Χ,Ĭ��Ϊ[0,1]
% v:SVM Cross Validation����,Ĭ��Ϊ5

%%
MAXGEN = gsa_option.maxgen;
NIND = gsa_option.sizepop;
NVAR = 3;  %��������������룬��������Ϊ60λȾɫ��
PRECI = 20;
IntervalCodeLength=size(intervals_index,1);
GGAP = gsa_option.ggap;
trace = zeros(MAXGEN,2);

FieldID = ...
[rep([PRECI],[1,NVAR]); ...
[gsa_option.cbound(1),gsa_option.gbound(1),gsa_option.pbound(1);gsa_option.cbound(2),gsa_option.gbound(2),gsa_option.pbound(2);];...
[0,0,0;0,0,0;0,1,1;1,1,1]];

%Create an initial population
Chrom = crtbp(NIND,IntervalCodeLength+NVAR*PRECI);%�볤Ϊ�������䳤��+60,Ⱦɫ��ǰ��Ϊ�������䣬����Ϊ60λ֧������������
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
    %��bestwaveindexs_temp
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
    %ObjVΪ����������Ӧ�Ⱥ���ֵ
    ObjV(nind,1) = svmtrain(train_label,train_data_best,cmd);
end
%������Ž⼰����ţ�YΪ���Ž⣬IΪȾɫ������Ⱥ�е����
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

%ģ���˻����ȷ��
INITTEMP=gsa_option.inittempFactor*(BadMSE-BestMSE);
%����ϵ��
LOWERTEMP=gsa_option.lowertempFactor;
% disp(INITTEMP);
% disp(LOWERTEMP);

steps=MAXGEN;
hwait=waitbar(0,'��ȴ�>>>>>>>>');%������
step=steps/100;
while 1
	%% ������
    if steps-gen<=1
        str1=['��', num2str(k),'��','�������'];
        waitbar(gen/steps,hwait,str1);
        pause(0.05);
    else
        PerStr=fix(gen/step);
        str2=['��', num2str(k),'��','���ڽ�����',num2str(PerStr),'%'];
        waitbar(gen/steps,hwait,str2);
        pause(0.05);
    end

    %������Ӧ��ֵ���㣬����¶Ȳ�������Ӧ�Ⱥ������
    for chromNum = 1:size(ObjV,1)
        FitnV(chromNum,1)=exp(-(ObjV(chromNum,1)-min(ObjV))/INITTEMP);    
    end
    %�Ŵ�����
    SelCh = select('rws',Chrom,FitnV,GGAP);
    SelCh = recombin('xovsp',SelCh,0.7);
    SelCh = mut(SelCh);
    %ģ���˻���������������ά����
    SelCh = simulatedAnneal_intervals_multibits(SelCh,INITTEMP,FieldID,train_label,train_data,v,intervals_index);  
    
 %���Ӵ�Ŀ�꺯��MSE���� 
    cg = bs2rv(SelCh(:,IntervalCodeLength+1:end),FieldID);
    for nind = 1:size(SelCh,1)
        cmd = ['-v ',num2str(v),' -c ',num2str(cg(nind,1)),' -g ',num2str(cg(nind,2)),' -p ',num2str(cg(nind,3)),' -s 3'];
        %��bestWaveIndexs_temp
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
    %���Ӵ�������һ����Ⱥ
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
  %����ʹsvm�Ĳ���c��С  
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
    %��ǰ�������Ӧ�Ⱥ�ƽ����Ӧ��
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
close(hwait);%�رս�����
