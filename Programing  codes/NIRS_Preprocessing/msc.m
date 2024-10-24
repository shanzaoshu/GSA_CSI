function [Xcorrect]=msc(X,Xref)
%  msc pretreate the samples X with the Multiplicative Scatter Correct
%  Input 
%         X:the matrix of the sample spectra to be Correct
%         Xref:the matrix of the sample spectra to be ref
%  Output
%         Xcorrect:the sample spectras was Corrected from the X  
%  Programmer: zhimin zhang @ central south university on dec 13 ,2007
%  Reference: Chemometrics and Intelligent Laboratory Systems 29 (1995) 233-241


[nRow nCol]=size(X);

vChannelWeights=ones(1,nCol);
nConditionNumber=10^12;
vRefSpectrum=mean(Xref);
DF=sum(vChannelWeights.^2); 


ZMod = []; 
ZMod = [ZMod ones(1,nCol)'];
ZMod = [ZMod vRefSpectrum'];
[nModRow,nModCol]=size(ZMod);
ZModW=ZMod.*(vChannelWeights'*ones(1,nModCol));
ZW=X.*(ones(nRow,1)*vChannelWeights);
ZZ=ZModW'*ZModW;
[u,S,v]=svd(ZZ);s=diag(S)'; 
sMinimum=s(1)/sqrt(nConditionNumber);
sCorrected=max(s,sMinimum ); 
ZZCorrected=u*(diag(sCorrected))*v';  
InvZZ = inv(ZZCorrected);
B = InvZZ*ZModW'*ZW';

ModelParam = B';

Xcorrect = X; 


p = ModelParam(:,1);
Xcorrect = Xcorrect-(p*ones(1,nModRow));

p = ModelParam(:,2);
Xcorrect = Xcorrect./(p*ones(nModRow,1)');


