 function [rmsecv] = pls_crossvalind_wavelength(train_data,train_label,wavelengthIndex,k,ncomp)

%rmsecv为交叉验证的均方根误差
%train_data：光谱数据
%train_label：为某一成分化学分析结果
%wavelengthIndex：特征波长编号
%k为Kfold中的参数k
%ncomp：默认主成分数
%% K-fold crossvalind
[m_size n]=size(train_data); %train_data为样本集合。每一行为一个观察样本
indices = crossvalind('Kfold',m_size,k); %产生k个fold，即indices里有等比例的1-k
rmse=[];
ncomp_temp=ncomp;
for i=1:k
    test=(indices==i); %逻辑判断，每次循环选取一个fold作为测试集
    train=~test; %取test的补集作为训练集，即剩下的k-1个fold
    data_train=train_data(train,:); %以上得到的数都为逻辑值，用于校正集样本数据的选取
    data_train_wavelength=data_train(:,wavelengthIndex);
    label_train=train_label(train,:); %label为样本类别标签，用于校正集样本标签的选取
    data_test=train_data(test,:); %选取测试集的样本数据
    data_test_wavelength=data_test(:,wavelengthIndex);
    label_test=train_label(test,:);%选取测试集的标签
    %% PLS模型
    %保证主成分数小于数据矩阵列数
    ncomp=ncomp_temp;
    if ncomp>size(data_train_wavelength,2)
        ncomp=size(data_train_wavelength,2);
    end
    if ncomp>0
        [m n_size]=size(data_train_wavelength);
        ab_train=[data_train_wavelength label_train];  %校正集样本数据和标签
        mu=mean(ab_train);sig=std(ab_train); %求均值和标准差
        rr=corrcoef(ab_train);   %求相关系数矩阵
        ab=zscore(ab_train); %数据标准化
        a=ab(:,[1:n_size]);b=ab(:,[n_size+1:end]);  %提出标准化后的自变量和因变量数据
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] =plsregress(a,b,ncomp);%主成分默认取10
        n=size(a,2); m=size(b,2);%n是自变量的个数,m是因变量的个数
        beta3(1,:)=mu(n+1:end)-mu(1:n)./sig(1:n)*BETA([2:end],:).*sig(n+1:end); %原始数据回归方程的常数项
        beta3([2:n+1],:)=(1./sig(1:n))'*sig(n+1:end).*BETA([2:end],:); %计算原始变量x1,...,xn的系数，每一列是一个回归方程R
        %     yhat=repmat(beta3(1,:),[size(a,1),1])+ab_train(:,[1:n])*beta3([2:end],:);  %求label_train的预测值
        ab_test=[data_test_wavelength label_test];%测试集样本数据和标签
        a1=ab_test(:,[1:n_size]);b1=ab_test(:,[n_size+1:end]); 
        yhat_test=repmat(beta3(1,:),[size(a1,1),1])+ab_test(:,[1:n])*beta3([2:end],:); %求label_test的预测值
        rmse=[rmse sqrt(sum((yhat_test-label_test).^2)/m_size*k)]; %求每次的rmsep
    end

end
   rmsecv=sum(rmse)/k;