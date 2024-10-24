function NewChrom = simulatedAnneal_intervals_multibits(OldChrom,temperature,FieldID,train_label,train_data,v,intervals_index)
      %% ģ���˻����
    [popSize,chromLength]=size(OldChrom);
    IntervalCodeLength=size(intervals_index,1);
    a=[];
    a=[a randperm(chromLength,ceil(IntervalCodeLength/10))]; %���������Ŷ�λ��Ϊ��Ӧ������볤�ȳ���10��ȡ��
    bitFlagIndex= a;
    NVAR = 3;%ǰ��Ϊ�������䣬����ΪSVM��3��������20λ�볤��
    PRECI = 20;
    %����⹹����ʽ
    bitFlag= 1;
    for metropolisNum= 1:popSize
         %��һλ����
        if bitFlag==0;
            mutBitNum=randi([1,chromLength],1,1);
            oldIndividual=OldChrom(metropolisNum,:);        
            newIndividual=oldIndividual;
            newIndividual(mutBitNum)=~newIndividual(mutBitNum);
        end
        %ÿ������λ����
        if bitFlag==1;
             mutBitNum2=randi([IntervalCodeLength+1,IntervalCodeLength+PRECI],1,1);
             mutBitNum3=randi([IntervalCodeLength+PRECI+1,IntervalCodeLength+2*PRECI],1,1);
             mutBitNum4=randi([IntervalCodeLength+2*PRECI+1,chromLength],1,1);

            oldIndividual=OldChrom(metropolisNum,:);        
            newIndividual=oldIndividual;
            newIndividual(bitFlagIndex)=~newIndividual(bitFlagIndex);%��λ�Ŷ�
            newIndividual(mutBitNum2)=~newIndividual(mutBitNum2);%һλ�Ŷ�
            newIndividual(mutBitNum3)=~newIndividual(mutBitNum3);
            newIndividual(mutBitNum4)=~newIndividual(mutBitNum4);
        end
        cg_old = bs2rv(oldIndividual(IntervalCodeLength+1:end),FieldID);
        cg_new = bs2rv(newIndividual(IntervalCodeLength+1:end),FieldID);
        
        cmd_old = ['-v ',num2str(v),' -c ',num2str(cg_old(1)),' -g ',num2str(cg_old(2)),' -p ',num2str(cg_old(3)),' -s 3'];
        %��bestWaveIndexs
        bestWaveIndexs_old=[];
        intervasStringOld=oldIndividual(1:IntervalCodeLength);
        for i=1:IntervalCodeLength
            if intervasStringOld(i)==1
                for j=intervals_index(i,1):intervals_index(i,2)
                   bestWaveIndexs_old=[bestWaveIndexs_old j];
                end
            end
        end
        train_data_best_old=train_data(:,bestWaveIndexs_old);   
        ObjV_old= svmtrain(train_label,train_data_best_old,cmd_old);
        
        cmd_new = ['-v ',num2str(v),' -c ',num2str(cg_new(1)),' -g ',num2str(cg_new(2)),' -p ',num2str(cg_new(3)),' -s 3'];
        %��bestWaveIndexs
        bestWaveIndexs_new=[];
        intervasStringNew=newIndividual(1:IntervalCodeLength);
        for i=1:IntervalCodeLength
            if intervasStringNew(i)==1
                for j=intervals_index(i,1):intervals_index(i,2)
                   bestWaveIndexs_new=[bestWaveIndexs_new j];
                end
            end
        end
        train_data_best_new=train_data(:,bestWaveIndexs_new);       
        ObjV_new= svmtrain(train_label,train_data_best_new,cmd_new);
        Chrom(metropolisNum,:)=oldIndividual;
        dt=ObjV_new-ObjV_old;
        if dt<0
           Chrom(metropolisNum,:)=newIndividual;
        else if rand<exp(-dt/temperature);
                 Chrom(metropolisNum,:)=newIndividual;
            end
        end    
         
end
    NewChrom=Chrom;
end