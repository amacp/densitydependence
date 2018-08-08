function []=DIModel()
warning ('off','all');
%Initial parameters
r=0.2;gamma=50;s=0.9;sigma=1;m=0.05;beta=1;epsilon=0.05;pBins=10;

nMax = ceil(1.5*(-1 + s + r*s)*gamma/(r*s));
tmax=500;Deltat=10;

TMtrx = MakeTransitionMatrix(nMax,r,gamma,s,sigma,epsilon,m);
iVec=initVec(nMax,2,2);
[Nmtrx,pmtrx,nAvgVec,pAvgVec]=BuildDynamicMtrx(nMax,pBins,TMtrx,iVec,tmax,Deltat);
figure(1)
clf;
hold on
imagesc(log(abs(Nmtrx)))
hold off
figure(2)
clf;
hold on
imagesc(log(abs(pmtrx)))
hold off
csvwrite('DIN.csv',abs(Nmtrx))
csvwrite('DIp.csv',abs(pmtrx))
csvwrite('DIavg.csv',[nAvgVec,pAvgVec])
end
function index=ele(nMax,i1, i2)
    index=i1*(nMax+1)+i2+1;
end
function TMtrx = MakeTransitionMatrix(nMax,r,gamma,s,sigma,epsilon,m)
%Birth matrix
bMatrix = zeros((nMax+1)^2,(nMax+1)^2);
for i1 =0:nMax
    for i2=0:nMax
        temp=0;
        for j1 = 0:nMax
            for j2=0:nMax
                if(j1>=i1 && j2>=i2)
                    if(i1*r*(1-(i1+i2)/gamma)>0) b1=exp(-i1*r*(1-(i1+i2)/gamma))*(i1*r*(1-(i1+i2)/gamma))^(j1-i1)/factorial(j1-i1);
                    elseif(j1==0) b1=1;
                    else b1=0;
                    end
                    if(i2*r*(1-(i1+i2)/gamma)>0) b2=exp(-i2*r*(1-(i1+i2)/gamma))*(i2*r*(1-(i1+i2)/gamma))^(j2-i2)/factorial(j2-i2);
                    elseif(j2==0) b2=1;
                    else b2=0;
                    end
                    bMatrix(ele(nMax,j1,j2),ele(nMax,i1,i2))=b1*b2;
                    temp=temp+bMatrix(ele(nMax,j1,j2),ele(nMax,i1,i2));
                end
            end
        end
        bMatrix(ele(nMax,i1,i2),ele(nMax,i1,i2))=bMatrix(ele(nMax,i1,i2),ele(nMax,i1,i2))+(1-temp);
    end
end
%Death matrix

dMatrix = zeros((nMax+1)^2,(nMax+1)^2);
for i1=0:nMax
    for i2=0:nMax
        for j1=0:nMax
            for j2=0:nMax
                if(j1<=i1 && j2<=i2)
                    dMatrix(ele(nMax,j1,j2),ele(nMax,i1,i2))=nchoosek(i1,j1)*s^j1*(1-s)^(i1-j1)*nchoosek(i2,j2)*(s-sigma*epsilon)^(j2)*(1-s+sigma*epsilon)^(i2-j2);
                end
            end
        end
    end
end

%Migration matrix

mMatrix = zeros((nMax+1)^2,(nMax+1)^2);
for i1=0:nMax
    for i2=0:nMax
        temp=0;
        for j1=0:nMax
            for j2=0:nMax
                if(j1>=i1 && j2>=i2)
                    mMatrix(ele(nMax,j1,j2),ele(nMax,i1,i2))=exp(-m)*m^(j1-i1)/factorial(j1-i1)*exp(-m)*m^(j2-i2)/factorial(j2-i2);
                    temp=temp+mMatrix(ele(nMax,j1,j2),ele(nMax,i1,i2));
                end
            end
        end
        mMatrix(ele(nMax,i1,i2),ele(nMax,i1,i2))=mMatrix(ele(nMax,i1,i2),ele(nMax,i1,i2))+(1-temp);
    end
end

TMtrx = mMatrix*dMatrix*bMatrix;
% figure(1)
% clf;
% hold on
% imagesc(flipud(log(abs(TMtrx))))
% hold off

end
function ivec=initVec(nMax,n01,n02)
    ivec = zeros((nMax+1)^2,1);
    ivec(ele(nMax,n01,n02),1)=1;
end
function [nTVec,pVec,nAvg,pAvg]=ToNTP(nMax,pBins,vecIn)
%This function takes in a vector in the form of TMtrx^t*iVec and outputs
%two vectors, a vector for the probabiltiy of nTot and a vector for the
%probabiltiy of p.
pVec=zeros(pBins,1);
nTVec=zeros(nMax+1,1);
nAvg=0;pAvg=0;
for i1=0:nMax
    for i2=0:nMax
      temp=i1/(i1+i2);
      nAvg=abs(nAvg+(i1+i2)*vecIn(ele(nMax,i1,i2)));
      if(i1+i2>0)
          pAvg=abs(pAvg+temp*vecIn(ele(nMax,i1,i2)));
      end
      for p=1:pBins
          if((p-1)/pBins<temp&&temp<=p/pBins)
              pVec(p)=pVec(p)+vecIn(ele(nMax,i1,i2));
          end
      end
      if(i1+i2<=nMax)
          nTVec(i1+i2+1)=nTVec(i1+i2+1)+vecIn(ele(nMax,i1,i2));
      else
          nTVec(nMax+1)=nTVec(nMax+1)+vecIn(ele(nMax,i1,i2));
      end
    end
end
end
function [Nmtrx,pmtrx,NAvgVec,pAvgVec]=BuildDynamicMtrx(nMax,pBins,TMtrx,iVec,tmax,Deltat)
    nSteps=floor(tmax/Deltat);
    [V,D]=eig(TMtrx);
    NAvgVec=zeros(nSteps+1,1);
    pAvgVec=zeros(nSteps+1,1);
    Nmtrx = zeros(nMax+1,nSteps+1);
    pmtrx = zeros(pBins,nSteps+1);
    [nTemp,pTemp,navg,pavg]=ToNTP(nMax,pBins,iVec);
    NAvgVec(1)=navg;pAvgVec(1)=pavg;
    Nmtrx(:,1)=nTemp;pmtrx(:,1)=pTemp;
    for t = 1:nSteps
        [nTemp,pTemp,navg,pavg]=ToNTP(nMax,pBins,V*(D.^(t*Deltat))/V*iVec);
        Nmtrx(:,t+1)=nTemp;pmtrx(:,t+1)=pTemp;
        NAvgVec(t+1)=navg;pAvgVec(t+1)=pavg;
    end
end