function []=EcologicalModel()
%Initial parameters
r=0.2;gamma=50;s=0.9;sigma=1;m=0.05;beta=1;epsilon=0.05;pBins=10;

Nmax = ceil(2*(((-1 + s + r*s)*gamma)/(r*s)));
tmax=200;

TMatrix = MakeTransitionMatrix(Nmax, r,gamma,s,m);
ivec=initVec(Nmax,5);
output = BuilddynamicMatrix(Nmax,TMatrix,ivec,tmax);
imagesc(flipud(log(abs(output))))
csvwrite('ecoN.csv',flipud(abs(output)))
end

function TMatrix = MakeTransitionMatrix(Nmax,r,gamma,s,m)

%Birth matrix
bMatrix = zeros(Nmax+1,Nmax+1);
for i = 0:Nmax
    temp = 0;
    lambda = i*r*(1-i/gamma);
    for j = 0:Nmax
        if(j>=i&&lambda>0)
            bMatrix(j+1,i+1) = (exp(-lambda)*lambda^(j-i))/factorial((j-i));
            temp = temp + bMatrix(j+1,i+1);
        end
    end
    bMatrix(i+1,i+1) = bMatrix(i+1,i+1) + (1-temp);
end
%Death matrix

dMatrix = zeros(Nmax+1,Nmax+1);

dMatrix(1,1)=1;
for i = 1:Nmax
    for j = 0:Nmax
        if(j<=i)
            dMatrix(j+1,i+1) = nchoosek(i,j)*s^j*(1-s)^(i-j);
        end
    end
end


%Migration matrix

mMatrix = zeros(Nmax+1,Nmax+1);

mMatrix(1,1)=1;
for i = 1:Nmax
    temp = 0;
    for j = 0:Nmax
        if(j>=i)
            mMatrix(j+1,i+1) = (exp(-m)*m^(j-i))/factorial((j-i));
            temp = temp + mMatrix(j+1,i+1);
        end
    end
    
    mMatrix(i+1,i+1) = mMatrix(i+1,i+1) + (1-temp);
end

TMatrix = mMatrix*dMatrix*bMatrix;

end

%Initial conditions

function ivec=initVec(Nmax,n0)
    ivec = zeros(Nmax+1,1);
    assert(n0 <= Nmax);
    assert(n0 >= 0);
    ivec(n0+1,1)=1;
end

function dynamicMatrix = BuilddynamicMatrix(nMax,TMatrix,ivec,tmax)
[V,D]=eig(TMatrix);
dynamicMatrix = zeros(nMax+1,tmax+1);
dynamicMatrix(:,1)=transpose(ivec);
for t = 1: tmax
    dynamicMatrix(:,t+1) = V*(D.^t)/V*ivec;
end
end

function expectedvec = Buildexpectedvector(dynamicMatrix,tmax,nMax)
temp = 0:nMax;
expectedvec = transpose(dynamicMatrix)*transpose(temp);
end