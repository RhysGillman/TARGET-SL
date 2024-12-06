function [ z_sub_network ] = sample_mutation_network( z,B_fda )
%Function:
%         forming the mutation network for each sample
%Input: 
%         z:original networks
%         mutation�������ı��������
%Output��
%         z_sub_network�������ı������������� 
mutation=B_fda;
N=max(max(z));
%�ڽӾ���
A=zeros(N);
for i=1:N
    A(z(i,2),z(i,1))=1;
end
%*************��������������mutation ����**************************
[ q_new ] = RWR( A,mutation );%������������㷨
%�����������Ľṹ�����¼���mutation gene�����������ߵĸ���
%��ͳ�������ԣ�ѡȡ����������Ľڵ㣬����sample-specific��mutation network
number=1;ind=1;
num_max=50;
while ind
dir_srand=dir_generate_srand(A);
[ q_rand ] = RWR( dir_srand,mutation );%������������㷨
Q_rand(:,number)=q_rand;
number=number+1;
if number>num_max
    ind=0;
end
end
%ͳ��100�β����ڵ������Ե�ƽ��ֵ�ͷ���
%������p_value
[row,colunm]=size(Q_rand);
for i=1:row
    p(i,1)=(max(Q_rand(i,:))<q_new(i,1));
end
p=double(p);
q_new0=q_new;
q_new0(q_new0>0.0001)=1;
p_value=q_new0.*p;
down_stream=find(p_value>0);
subnetwork_nodes=unique([mutation;down_stream]');
%��ȡ��������Ľṹ
[x1,y1]=ismember(z(:,1),subnetwork_nodes);
[x2,y2]=ismember(z(:,2),subnetwork_nodes);
ind=double(x1).*double(x2);
z_sub_network=z(find(ind~=0),:);


end

