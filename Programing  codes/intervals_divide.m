function intervals_index = intervals_divide(nirs_data,intervals)
%nirs_dataΪ�������ݼ�
%intervalsΪ�����������
%intervals_indexΪ�����ţ�intervals�У�2�У���ʼ��źͽ�����ţ�
[nint,mint]=size(intervals);
[n,m]=size(nirs_data);
if mint > 1
    allint=[(1:round(mint/2)+1)' [intervals(1:2:mint)';1] [intervals(2:2:mint)';m]];
    intervals=round(mint/2);
    intervalsequi=0;
else
    vars_left_over=mod(m,intervals);
    N=fix(m/intervals);
    % Distributes vars_left_over in the first "vars_left_over" intervals
    startint=[(1:(N+1):(vars_left_over-1)*(N+1)+1)'; ((vars_left_over-1)*(N+1)+1+1+N:N:m)'];
    endint=[startint(2:intervals)-1; m];
    allint=[(1:intervals+1)' [startint;1] [endint;m]];
end
intervals_index=allint(1:intervals,2:3);