function [BestWavelengthIndexs,BestIntervalsString,BestMSE] = gsaPLSintervalsForRegress(train_data,train_label,gsa_option,intervals_index,k)
%% ������ʼ��
%train_labelΪĳһ��һ�ɷ�ֵ���Ը�ֵΪ��׼�����Ż�
if nargin == 2
    gsa_option = struct('inittempFactor',100,'lowertempFactor',0.8,'maxgen',200,'sizepop',30,'ggap',0.9,'m',10,'v',5);
end
% inittemp:��ʼ�¶ȣ�Ĭ��Ϊ100
%lowertemp:����ϵ����Ĭ��Ϊ0.8
% maxgen:���Ľ�������,Ĭ��Ϊ200,һ��ȡֵ��ΧΪ[100,500]
% sizepop:��Ⱥ�������,Ĭ��Ϊ20,һ��ȡֵ��ΧΪ[20,100]
% v:SVM Cross Validation����,Ĭ��Ϊ5

MAXGEN = gsa_option.maxgen;
NIND = gsa_option.sizepop;
IntervalCodeLength=size(intervals_index,1);
GGAP = gsa_option.ggap;
trace = zeros(MAXGEN,4);%���α������Ŀ�꺯����ƽ��Ŀ�꺯���������Ӧ�Ⱥ�������С��Ӧ�Ⱥ���ֵ
v = gsa_option.v;
m =gsa_option.m;

%% Create an initial population
Chrom = crtbp(NIND,IntervalCodeLength);


gen = 1;
BestMSE = inf; %inf��ʾ�����
BadMSE =inf;

%% ������ȺĿ�꺯��ֵ
for nind = 1:NIND
    %ObjVΪ����������Ӧ�Ⱥ���ֵRMSECV,Ҫ����ÿһ��Ⱦɫ���RMSE
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
%������Ž⼰����ţ�YΪ���Ž⣬IΪȾɫ������Ⱥ�е����
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
%ģ���˻����ȷ��
INITTEMP=gsa_option.inittempFactor*(BadMSE-BestMSE);
%����ϵ��
LOWERTEMP=gsa_option.lowertempFactor;
%% 
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
  
    %����¶Ȳ�������Ӧ�Ⱥ������
    for chromNum = 1:size(ObjV,1)
        FitnV(chromNum,1)=exp(-(ObjV(chromNum,1)-min(ObjV))/INITTEMP);    
    end
    %�Ŵ�����
    SelCh = select('rws',Chrom,FitnV,GGAP);%����ѡ��
   % SelCh = recombin('xovsp',SelCh,0.7); %���㽻��
    SelCh = recombin('recdis',SelCh,0.7); %��ɢ����
    SelCh = mut(SelCh,0.01);%��ɢ����
    %ģ���˻����
    SelCh = pls_interval_simulatedAnneal(SelCh,INITTEMP,train_data,train_label,m,intervals_index,v);  
    
 %���Ӵ�Ŀ�꺯��MSE���� 
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
    %���Ӵ�������һ����Ⱥ
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
      %��ǰ�������Ӧ�Ⱥ�ƽ����Ӧ��
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

  close(hwait);%�رս�����
%%
% figure;
% hold on;
% set(gcf,'color','white');
% 
% %���û�ͼ��С
% set(gcf,'Units','centimeters','Position',[10 10 17 13]);
% 
% trace = round(trace*10000)/10000;
% plot(trace(1:gen,1),'k*-');
% plot(trace(1:gen,2),'kd-');
% 
% tl=legend('���Ŀ�꺯��ֵ','ƽ��Ŀ�꺯��ֵ');
% %����tl�������С
% set(tl,'fontsize',9);
% xlabel('��������','FontSize',9);
% ylabel('Ŀ�꺯��ֵ','FontSize',9);
% %����������̶�ֵ�����С
% 
% figure;
% hold on;
% set(gcf,'color','white');
% 
% %���û�ͼ��С
% set(gcf,'Units','centimeters','Position',[10 10 17 13]);
% 
% trace = round(trace*10000)/10000;
% plot(trace(1:gen,3),'k*-');
% plot(trace(1:gen,4),'kd-');
% 
% tl=legend('ƽ����Ӧ��','��С��Ӧ��');
% %����tl�������С
% set(tl,'fontsize',9);
% xlabel('��������','FontSize',9);
% ylabel('��Ӧ�Ⱥ���ֵ','FontSize',9);
% %����������̶�ֵ�����С












