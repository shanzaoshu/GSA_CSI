 function [rmsecv] = pls_crossvalind_wavelength(train_data,train_label,wavelengthIndex,k,ncomp)

%rmsecvΪ������֤�ľ��������
%train_data����������
%train_label��Ϊĳһ�ɷֻ�ѧ�������
%wavelengthIndex�������������
%kΪKfold�еĲ���k
%ncomp��Ĭ�����ɷ���
%% K-fold crossvalind
[m_size n]=size(train_data); %train_dataΪ�������ϡ�ÿһ��Ϊһ���۲�����
indices = crossvalind('Kfold',m_size,k); %����k��fold����indices���еȱ�����1-k
rmse=[];
ncomp_temp=ncomp;
for i=1:k
    test=(indices==i); %�߼��жϣ�ÿ��ѭ��ѡȡһ��fold��Ϊ���Լ�
    train=~test; %ȡtest�Ĳ�����Ϊѵ��������ʣ�µ�k-1��fold
    data_train=train_data(train,:); %���ϵõ�������Ϊ�߼�ֵ������У�����������ݵ�ѡȡ
    data_train_wavelength=data_train(:,wavelengthIndex);
    label_train=train_label(train,:); %labelΪ��������ǩ������У����������ǩ��ѡȡ
    data_test=train_data(test,:); %ѡȡ���Լ�����������
    data_test_wavelength=data_test(:,wavelengthIndex);
    label_test=train_label(test,:);%ѡȡ���Լ��ı�ǩ
    %% PLSģ��
    %��֤���ɷ���С�����ݾ�������
    ncomp=ncomp_temp;
    if ncomp>size(data_train_wavelength,2)
        ncomp=size(data_train_wavelength,2);
    end
    if ncomp>0
        [m n_size]=size(data_train_wavelength);
        ab_train=[data_train_wavelength label_train];  %У�����������ݺͱ�ǩ
        mu=mean(ab_train);sig=std(ab_train); %���ֵ�ͱ�׼��
        rr=corrcoef(ab_train);   %�����ϵ������
        ab=zscore(ab_train); %���ݱ�׼��
        a=ab(:,[1:n_size]);b=ab(:,[n_size+1:end]);  %�����׼������Ա��������������
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] =plsregress(a,b,ncomp);%���ɷ�Ĭ��ȡ10
        n=size(a,2); m=size(b,2);%n���Ա����ĸ���,m��������ĸ���
        beta3(1,:)=mu(n+1:end)-mu(1:n)./sig(1:n)*BETA([2:end],:).*sig(n+1:end); %ԭʼ���ݻع鷽�̵ĳ�����
        beta3([2:n+1],:)=(1./sig(1:n))'*sig(n+1:end).*BETA([2:end],:); %����ԭʼ����x1,...,xn��ϵ����ÿһ����һ���ع鷽��R
        %     yhat=repmat(beta3(1,:),[size(a,1),1])+ab_train(:,[1:n])*beta3([2:end],:);  %��label_train��Ԥ��ֵ
        ab_test=[data_test_wavelength label_test];%���Լ��������ݺͱ�ǩ
        a1=ab_test(:,[1:n_size]);b1=ab_test(:,[n_size+1:end]); 
        yhat_test=repmat(beta3(1,:),[size(a1,1),1])+ab_test(:,[1:n])*beta3([2:end],:); %��label_test��Ԥ��ֵ
        rmse=[rmse sqrt(sum((yhat_test-label_test).^2)/m_size*k)]; %��ÿ�ε�rmsep
    end

end
   rmsecv=sum(rmse)/k;