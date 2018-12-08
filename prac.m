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
% % %          data=sum(sum(number([1 5 6 7 8 9 10 12 17],4:12)));  %一个月里的高科技岗位
% %          Data_all=[Data_all;data];%36个月的高科技岗位
% %           mdata=sum(sum(number(:,4:end)));
% %           Data=[Data,mdata];%计算每个月岗位需求量
% 
%     end
% end
% save Data_all
% save Data


% talent=sum(Data_all(:,2:end),2);
% low=Data_all(:,[2 3 4]);mid=Data_all(:,[5,6]);high=Data_all(:,[7 8 9]);%对学历进行分类
% level=sum(low,2).*1+sum(mid,2).*2+sum(high,2).*3;%加权 1，2,3
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

% %%%%%%%%%%%%%%%%%%%%%第一题模型求解
% X=[ones(size(talent)) level job expect ];
% [b,bint,r,rint,stats] = regress(talent,X)
% figure(4) %拟合后设置属性
%  title('The relationship between educational background and talent demand')
% figure(5)
%  xlabel('Job demand'),ylabel('Talent demand')
%  title('The relationship between Job demand and talent demand')
%  figure(1)
%  title('The relationship between Desire profession and talent demand')
%  xlabel('Desire profession'),ylabel('Talent demand')
% 
% 
% %%%%%%%%%%%%%%%%%%%%%第二题
% data_1=sum(Data(1:12));
% data_2=sum(Data(13:24));
% data_3=sum(Data(25:36));
% 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ARIMA
s = 12; %周期是24
x = Data';%初始数据的录入
n = 36; %预报的个数
m1 = length(x); %原始的数据的个数
for i = s+1:m1
    y(i-s) = x(i) - x(i-s);%进行周期差分变换
end

w = diff(y); %消除趋势性的差分运算
m2 = length(w);k=0;
for i = 0:2
    for j = 0:2
        if i == 0 & j == 0
            continue
        elseif i == 0
            ToEstMd = arima('MALags',1:j,'Constant',0); %指定模型的结构
        elseif j == 0
            ToEstMd = arima('ARLags',1:i,'Constant',0); %指定模型的结构
        else
            ToEstMd = arima('ARLags',1:i,'MALags',1:j,'Constant',0); %指定模型的结构
        end
        k = k + 1;
        R(k) = i;
        M(k) = j;
        [EstMd,EstParamCov,LogL,info] = estimate(ToEstMd,w');%模型拟合
        numParams = sum(any(EstParamCov));%计算拟合参数的个数
        [aic(k),bic(k)] = aicbic(LogL,numParams,m2);
    end
end
fprintf('R,M,AIC,BIC的对应值如下\n%f');%显示计算结果
check  = [R',M',aic',bic']


r=input('输入阶数R=');m=input('输入阶数M=');
ToEstMd = arima('ARLags',1:r,'MALags',1:m,'Constant',0);%指定模型的结构
[EstMd,EstParamCov,LogL,info] = estimate(ToEstMd,w');%模型拟合 
w_Forecast = forecast(EstMd,n,'Y0',w');
yhat = y(end) + cumsum(w_Forecast); %一阶差分的还原值
for j = 1:n
    x(m1 + j) = yhat(j) + x(m1+j-s); %x的预测值
end
pre=x(m1+1:end)
new=[x;pre];
newy=[sum(new(1:12)) sum(new(13:24)) sum(new(25:36)) sum(new(37:48)) sum(new(49:60)) sum(new(61:72))];



% %%%%%%%%%%%%问题二绘图
% 
% bar(Data,'r')
% xlabel('Month'),ylabel('Posts')  %36个月的柱状图
% subplot(1,2,1)
% autocorr(x)
% subplot(1,3,3)
% parcorr(x)
% 
% bar(y,'r')
% xlabel('Month'),ylabel('Posts')  %36个月的柱状图
% subplot(1,2,1)
% autocorr(y)
% subplot(1,2,2)
% parcorr(y)
% 
% bar(w,'r')
% xlabel('Month'),ylabel('Posts')  %36个月的柱状图
% subplot(1,2,1)
% autocorr(w)
% subplot(1,2,2)
% parcorr(w)
% bar([x;ans])
% xlabel('Momth');ylabel('Posts')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%问题四
% d2014=7270000*(1-0.792-0.014-0.037-0.026);
% d2015=7490000*(1-0.774-0.015-0.039-0.026);
% d2016=7560000*(1-0.773-0.015-0.04-0.024);
% d2017=7950000*(1-0.771-0.016-0.034-0.024);
% syms a b;
% c = [a b]';
% 
% %原始数列 A
% A = [d2014 d2015 d2016 d2017];
% n = length(A);
% 
% %对原始数列 A 做累加得到数列 B
% B = cumsum(A);
% 
% %对数列 B 做紧邻均值生成
% for i = 2:n
%     C(i) = (B(i) + B(i - 1))/2; 
% end
% C(1) = [];
% 
% %构造数据矩阵 
% B = [-C;ones(1,n-1)];
% Y = A; Y(1) = []; Y = Y';
% 
% %使用最小二乘法计算参数 a(发展系数)和b(灰作用量)
% c = inv(B*B')*B*Y;
% c = c';
% a = c(1); b = c(2);
% 
% %预测后续数据
% F = []; F(1) = A(1);
% for i = 2:(n+10)
%     F(i) = (A(1)-b/a)/exp(a*(i-1))+ b/a;
% end
% 
% %对数列 F 累减还原,得到预测出的数据
% G = []; G(1) = A(1);
% for i = 2:(n+10)
%     G(i) = F(i) - F(i-1); %得到预测出来的数据
% end
% 
% disp('预测数据为：');
% G
% 
% %模型检验
% 
% H = G(1:4);
% %计算残差序列
% epsilon = A - H;
% 
% %法一：相对残差Q检验
% %计算相对误差序列
% delta = abs(epsilon./A);
% %计算相对误差Q
% disp('相对残差Q检验：')
% Q = mean(delta)
% 
% %法二：方差比C检验
% disp('方差比C检验：')
% C = std(epsilon, 1)/std(A, 1)
% 
% %法三：小误差概率P检验
% S1 = std(A, 1);
% tmp = find(abs(epsilon - mean(epsilon))< 0.6745 * S1);
% disp('小误差概率P检验：')
% P = length(tmp)/n
% % 
% % %绘制曲线图
% % t1 = 2014:2017;
% % t2 = 2014:2027;
% % figure(2)
% % plot(t1, A,'ro'); hold on;
% % plot(t2, G, 'g-');
% % % xlabel('年份'); ylabel('污水量/亿吨');
% % % legend('实际污水排放量','预测污水排放量');
% % % title('长江污水排放量增长曲线');
% % grid on;
% % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%问题三

% [data1,data2,data3,data4,data5]=textread('C:\Users\DELL\Desktop\APMCM代码\alldata.txt','%n%n%n%n%n',30);
% x=[data1./100000,data2./100000,data3./100000,data4./100000];
% x(:,1)=[]
% 





