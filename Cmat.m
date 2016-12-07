load K_M_Newmark.mat
  OPTS.issym=1;
  OPTS.isreal=1;
[minvals,minvallocs]=sort(diag(Kr)./diag(Mr));
  
  shift=minvals(min(7,length(minvals)));
  [fms,f]=eigs((Kr+Kr')/2+shift*(Mr+Mr')/2,(Mr+Mr')/2,min([ size(Kr,1) ...
		    max([floor(sqrt(size(Kr,1))) 250])]),0,OPTS);
  fs=sqrt(diag(f)-shift)/2/pi;

  [fs,fsindex]=sort(fs);
  fms=fms(:,fsindex);

for i=1:150    
% This loop is for a diagonal matrix 2*Xi*wi
    ZW(i,i)=2*0.02*fs(i)*2*pi; 
end

C=(fms')\ZW/fms
save('C.mat','C')
