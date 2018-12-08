% clc,clear
% path = ('C:\Annex Job market of A-City Market Demand Statistics');
% file1 = dir(path);
% Data_all=[];
% Data=[];
% for i = 4:length(file1)
%   name1 = strcat(path,'\',file1(i).name); 
%     file2=dir(name1);
%   
%     for j=3:length(file2)
%         name2=strcat(name1,'\',file2(j).name);
%       
%          [number,~,~]=xlsread(name2); 
%          data=number(:,3:end);
%           Data_all=[Data_all;data];
% % %          if (i==5&j>5)||(i>5)
% % %          number(1,:)=[];
% % %         
% % %          end
% % %          data=sum(sum(number([1 5 6 7 8 9 10 12 17],4:12)));  %һ������ĸ߿Ƽ���λ
% %          Data_all=[Data_all;data];%36���µĸ߿Ƽ���λ
% %           mdata=sum(sum(number(:,4:end)));
% %           Data=[Data,mdata];%����ÿ���¸�λ������
% 
%     end
% end
% save Data_all
% save Data


% talent=sum(Data_all(:,2:end),2);
% low=Data_all(:,[2 3 4]);mid=Data_all(:,[5,6]);high=Data_all(:,[7 8 9]);%��ѧ�����з���
% level=sum(low,2).*1+sum(mid,2).*2+sum(high,2).*3;%��Ȩ 1��2,3
% figure (1)
% plot(talent,level,'.')
% xlabel('Desire educaional background'),ylabel('Talent demand')
% % title('The relationship between educational background and talent demand')
% job=Data_all(:,end);
% figure (2)
% plot(talent,job,'.k')
% xlabel('Job demand'),ylabel('Talent demand')
% % title('The relationship between Job demand and talent demand')
% expect=Data_all(:,1);
% figure (3)
% plot(talent,expect,'.r')
% % set(gca,'XTick',[0:100:1000]) 
% xlabel('Desire profession'),ylabel('Talent demand')
% title('The relationship between Desire profession and talent demand')

% %%%%%%%%%%%%%%%%%%%%%��һ��ģ�����
% X=[ones(size(talent)) level job expect ];
% [b,bint,r,rint,stats] = regress(talent,X)
% figure(4) %��Ϻ���������
%  title('The relationship between educational background and talent demand')
% figure(5)
%  xlabel('Job demand'),ylabel('Talent demand')
%  title('The relationship between Job demand and talent demand')
%  figure(1)
%  title('The relationship between Desire profession and talent demand')
%  xlabel('Desire profession'),ylabel('Talent demand')
% 
% 
% %%%%%%%%%%%%%%%%%%%%%�ڶ���
% data_1=sum(Data(1:12));
% data_2=sum(Data(13:24));
% data_3=sum(Data(25:36));
% 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ARIMA
s = 12; %������24
x = Data';%��ʼ���ݵ�¼��
n = 36; %Ԥ���ĸ���
m1 = length(x); %ԭʼ�����ݵĸ���
for i = s+1:m1
    y(i-s) = x(i) - x(i-s);%�������ڲ�ֱ任
end

w = diff(y); %���������ԵĲ������
m2 = length(w);k=0;
for i = 0:2
    for j = 0:2
        if i == 0 & j == 0
            continue
        elseif i == 0
            ToEstMd = arima('MALags',1:j,'Constant',0); %ָ��ģ�͵Ľṹ
        elseif j == 0
            ToEstMd = arima('ARLags',1:i,'Constant',0); %ָ��ģ�͵Ľṹ
        else
            ToEstMd = arima('ARLags',1:i,'MALags',1:j,'Constant',0); %ָ��ģ�͵Ľṹ
        end
        k = k + 1;
        R(k) = i;
        M(k) = j;
        [EstMd,EstParamCov,LogL,info] = estimate(ToEstMd,w');%ģ�����
        numParams = sum(any(EstParamCov));%������ϲ����ĸ���
        [aic(k),bic(k)] = aicbic(LogL,numParams,m2);
    end
end
fprintf('R,M,AIC,BIC�Ķ�Ӧֵ����\n%f');%��ʾ������
check  = [R',M',aic',bic']


r=input('�������R=');m=input('�������M=');
ToEstMd = arima('ARLags',1:r,'MALags',1:m,'Constant',0);%ָ��ģ�͵Ľṹ
[EstMd,EstParamCov,LogL,info] = estimate(ToEstMd,w');%ģ����� 
w_Forecast = forecast(EstMd,n,'Y0',w');
yhat = y(end) + cumsum(w_Forecast); %һ�ײ�ֵĻ�ԭֵ
for j = 1:n
    x(m1 + j) = yhat(j) + x(m1+j-s); %x��Ԥ��ֵ
end
pre=x(m1+1:end)
new=[x;pre];
newy=[sum(new(1:12)) sum(new(13:24)) sum(new(25:36)) sum(new(37:48)) sum(new(49:60)) sum(new(61:72))];



% %%%%%%%%%%%%�������ͼ
% 
% bar(Data,'r')
% xlabel('Month'),ylabel('Posts')  %36���µ���״ͼ
% subplot(1,2,1)
% autocorr(x)
% subplot(1,3,3)
% parcorr(x)
% 
% bar(y,'r')
% xlabel('Month'),ylabel('Posts')  %36���µ���״ͼ
% subplot(1,2,1)
% autocorr(y)
% subplot(1,2,2)
% parcorr(y)
% 
% bar(w,'r')
% xlabel('Month'),ylabel('Posts')  %36���µ���״ͼ
% subplot(1,2,1)
% autocorr(w)
% subplot(1,2,2)
% parcorr(w)
% bar([x;ans])
% xlabel('Momth');ylabel('Posts')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%������
% d2014=7270000*(1-0.792-0.014-0.037-0.026);
% d2015=7490000*(1-0.774-0.015-0.039-0.026);
% d2016=7560000*(1-0.773-0.015-0.04-0.024);
% d2017=7950000*(1-0.771-0.016-0.034-0.024);
% syms a b;
% c = [a b]';
% 
% %ԭʼ���� A
% A = [d2014 d2015 d2016 d2017];
% n = length(A);
% 
% %��ԭʼ���� A ���ۼӵõ����� B
% B = cumsum(A);
% 
% %������ B �����ھ�ֵ����
% for i = 2:n
%     C(i) = (B(i) + B(i - 1))/2; 
% end
% C(1) = [];
% 
% %�������ݾ��� 
% B = [-C;ones(1,n-1)];
% Y = A; Y(1) = []; Y = Y';
% 
% %ʹ����С���˷�������� a(��չϵ��)��b(��������)
% c = inv(B*B')*B*Y;
% c = c';
% a = c(1); b = c(2);
% 
% %Ԥ���������
% F = []; F(1) = A(1);
% for i = 2:(n+10)
%     F(i) = (A(1)-b/a)/exp(a*(i-1))+ b/a;
% end
% 
% %������ F �ۼ���ԭ,�õ�Ԥ���������
% G = []; G(1) = A(1);
% for i = 2:(n+10)
%     G(i) = F(i) - F(i-1); %�õ�Ԥ�����������
% end
% 
% disp('Ԥ������Ϊ��');
% G
% 
% %ģ�ͼ���
% 
% H = G(1:4);
% %����в�����
% epsilon = A - H;
% 
% %��һ����Բв�Q����
% %��������������
% delta = abs(epsilon./A);
% %����������Q
% disp('��Բв�Q���飺')
% Q = mean(delta)
% 
% %�����������C����
% disp('�����C���飺')
% C = std(epsilon, 1)/std(A, 1)
% 
% %������С������P����
% S1 = std(A, 1);
% tmp = find(abs(epsilon - mean(epsilon))< 0.6745 * S1);
% disp('С������P���飺')
% P = length(tmp)/n
% % 
% % %��������ͼ
% % t1 = 2014:2017;
% % t2 = 2014:2027;
% % figure(2)
% % plot(t1, A,'ro'); hold on;
% % plot(t2, G, 'g-');
% % % xlabel('���'); ylabel('��ˮ��/�ڶ�');
% % % legend('ʵ����ˮ�ŷ���','Ԥ����ˮ�ŷ���');
% % % title('������ˮ�ŷ�����������');
% % grid on;
% % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%������

% [data1,data2,data3,data4,data5]=textread('C:\Users\DELL\Desktop\APMCM����\alldata.txt','%n%n%n%n%n',30);
% x=[data1./100000,data2./100000,data3./100000,data4./100000];
% x(:,1)=[]
% 





