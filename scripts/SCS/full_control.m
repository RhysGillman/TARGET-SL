function [ manner,Nd ] = full_control( z )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
Z=[];
N=max(max(z));
Z(:,1)=z(:,1);Z(:,2)=z(:,2)+N;
all_B=[1:N]';
B=[];
[c,d]=size(Z);
for i=1:c
    B(Z(i,1),Z(i,2))=1;
    B(Z(i,2),Z(i,1))=1;
end
B=sparse(B);
% ��Ŀ��ڵ��йص�ƥ��
F = matching(B);
%����Markovʶ��ͬ��ƥ��
%MΪ����matchingʶ��ĳ�ʼƥ�䣬resultΪ����markovʶ��Ĳ�ͬƥ��
f=[];
f(:,1)=[1:N]';
f(:,2)=F(1:N,1);
PP=[];kk=1;
[row,colunm]=size(f);
for i=1:row
    if f(i,2)~=0
        PP(kk,1)=f(i,1);PP(kk,2)=f(i,2)-N;
        kk=kk+1;
    end
end
M=zeros(N);
[oo1,oo2]=size(PP);
for i=1:oo1
    M(PP(i,1),PP(i,2))=1;
end
manner=PP;
Nd=setdiff(all_B,PP(:,2));


end

