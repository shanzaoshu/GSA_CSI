%%�����ѡȡһ���������������ѵ����������ѡȡ��ȫ�������û�й��ɵĻ�����ѭ�򵥵Ĺ���
%%ѵ������ɷ����򵥣�����Ҫ����������ѡ����ÿ�����ѵ�������������ܲ���ܴ󣬲��ܱ�֤��ѡ���������Լ�ģ�͵�����������
function m = rs(X,N)

% Random Sampling Algorithm for selection of samples
% m = rs(X,N);
%
% X --> Matrix of instrumental responses
% N --> Number of samples to be selected 
%
% m --> Indexes of the selected samples

M = size(X,1); % Number of rows in X (samples)

rand_samples = randperm(M);

m = rand_samples(1:N);
