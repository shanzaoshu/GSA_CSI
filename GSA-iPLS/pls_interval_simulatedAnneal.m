function NewChrom = pls_interval_simulatedAnneal(OldChrom,temperature,train_data,train_label,m_bit,intervals_index,v)
%m_bit为扰动解的变异位数      
%% 模拟退火操作
    [popSize chromLength]=size(OldChrom);
    IntervalCodeLength=size(intervals_index,1);
    a=[];
    a=[a randperm(chromLength,m_bit)];
    %邻域解构建方式
    bitFlagIndex= a;
    for metropolisNum= 1:popSize
        oldIndividual=OldChrom(metropolisNum,:);        
        newIndividual=oldIndividual;
        newIndividual(bitFlagIndex)=~newIndividual(bitFlagIndex); %生成染色体的扰动解，对应bitFlagIndex中的基因位取反
        oldWavelengthIndexs=[];
        newWavelengthIndexs=[];
        %提取特征谱区对应的特征波长变量
        for i=1:IntervalCodeLength
            if oldIndividual(i)==1
                for j=intervals_index(i,1):intervals_index(i,2)
                    oldWavelengthIndexs=[oldWavelengthIndexs j];
                end
            end
            if newIndividual(i)==1
                for j=intervals_index(i,1):intervals_index(i,2)
                    newWavelengthIndexs=[newWavelengthIndexs j];
                end
            end
        end 
        ObjV_old=pls_crossvalind_wavelength(train_data,train_label,oldWavelengthIndexs,v,10);
        ObjV_new=pls_crossvalind_wavelength(train_data,train_label,newWavelengthIndexs,v,10);
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
