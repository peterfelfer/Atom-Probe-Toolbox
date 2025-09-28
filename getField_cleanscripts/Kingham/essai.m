clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%Model  Variable%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_resolution=1;%field resolution as 1 V/nm
Fmax=200;      %field maximum of model as 200 V/nm
z0_resolution=10000000000; %field distance resolution as 10000000000au 
zmax=1000000000000; %field maximum distance of model as 1000000000000au
zmin=1;
Lamda=0.8;     %screen length found by Lang and Kohn (p.279)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Import Table of Elemental constant%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IonizationTable=importdata('IonizationTable.mat'); %import elemental ionisation energy table
ConstantTable=importdata('ConstantTable.mat'); %import rest elemental variables such as mass, work energy and quantum number 
ElementList=importdata('ElementList.mat');  %import elemental name list
SubNumber=1;
figure('Position',[1,1,1440,900])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for e=11;%:1:3
%     size(ElementList,1)
                %%%%%%%Assign varables from IonizationTable%%%%%%%%%%%%%%%%%%
                I=zeros(1,1);
                I_raw=IonizationTable(e,:);
                i=2;
                while (isnan(I_raw(1,i-1))==0)
                   I(1,i)=I_raw(1,i-1);
                   i=i+1;
                end
                BB=I;
                BB=BB/27.2;               
                I=I/27.2;
                %%%%%%%Assign varables from ConstantTable%%%%%%%%%%%%%%%%%%
                m=ConstantTable(e,2);
                mm=ConstantTable(e,1);
                Phi=ConstantTable(e,3)/27.2;
                ElementName=ElementList(e,1);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                z0=(ones(1,Fmax/F_resolution)'*(zmin:z0_resolution:zmax))';
                z0=z0;
                zsize=size(z0,1);
                F=ones(1,zsize)'*(F_resolution:F_resolution:Fmax);
                F=F/514.2;
                Fsize=Fmax/F_resolution;
                E_stark=zeros(zsize,Fsize);
                E_image=zeros(zsize,Fsize);
                E_stark_2=zeros(zsize,Fsize);
                E_image_2=zeros(zsize,Fsize);
                D_E_stark=zeros(zsize,Fsize);
                D_E_image=zeros(zsize,Fsize);
                Zc=zeros(zsize,Fsize);
                Zz0=zeros(zsize,Fsize);

                SummationOfE_stark=zeros(zsize,Fsize,7);
                SummationOfE_image=zeros(zsize,Fsize,7);
                SummationOfZc=zeros(zsize,Fsize,7);
                SummationOfZz0=zeros(zsize,Fsize,7);
            for n=1:1:(size(I,2)-1)
                for f=1:1:(Fmax)
                    for z=1:1:zsize                    
                            E_stark(z,f)=(-1).*(17*m^2+10)/64*(F(z,f)/I(n+1))^2;    
                            E_image(z,f)=0.5*(n+0.5)/(z0(z,f)+Lamda);
                            if n==0
                                E_stark_2(z,f)=0; 
                                E_image_2(z,f)=0;
                            else
                                E_stark_2(z,f)=(-1).*(17*m^2+10)/64*(F(z,f)/I(n))^2; 
                                E_image_2(z,f)=0.5*(n+0.5-1)/(z0(z,f)+Lamda);
                            end
                            D_E_stark(z,f)=E_stark_2(z,f)-E_stark(z,f);
                            D_E_image(z,f)=E_stark_2(z,f)-E_image(z,f);            
                            Zc(z,f)=(I(n+1)-E_stark(z,f)-E_image(z,f)-Phi-Lamda*F(z,f))/F(z,f);
                            Zz0(z,f)=n+1+1+4.5/z0(z,f);
                    end 
                end
                SummationOfZc(:,:,n)=Zc(:,:);
                SummationOfZz0(:,:,n)=Zz0(:,:);
            n=n
            end
            Rest=zeros(zsize,Fsize);
            u=zeros(zsize,Fsize);
            R=zeros(zsize,Fsize);
            Ru=zeros(zsize,Fsize);
            term=zeros(zsize,Fsize);
            Pt=zeros(1,Fsize);
            SummationPt=zeros(1,Fsize);
            for n=1:1:(size(I,2)-1)
                A2v=I(n+1)/(6*pi*m*exp(2/3));
                B=BB(n+1);
                for f=1:1:Fsize
                    for z=1:1:zsize
                        Summation1=0; 
                        Summation2=0; 
                        Rest(z,f)=((n))*F(z,f).*(z0(z,f)+Lamda)+((n)^2)/(4*(z0(z,f)+Lamda))-(F(z,f))^(0.5);
                        for nn=1:1:(n)  
                            if nn==1
                               Zc2=0;
                               Summation1=0; 
                               Summation2=0;
                            else
                               Zc2=SummationOfZc(z,f,nn+1);
                               Summation1=Summation1+F(z,f).*(Zc2-SummationOfZc(z,f,nn)+Lamda); 
                               Summation2=Summation2+(2*(nn)+1)/(4*(Zc2-SummationOfZc(z,f,nn)+Lamda));
                            end  
                        end
                        u(z,f)=((Rest(z,f)-Summation1-Summation2)/mm*2)^0.5;
                        term(z,f)=(B-SummationOfZz0(z,f,n)*F(z,f)/B-F(z,f)*z0(z,f));
                        if term(z,f)<0
                            term(z,f)=0;
                        else
                            ;
                        end
                        R(z,f)=6*pi*A2v*F(z,f)/((2^(5/2)).*(B^(3/2)-term(z,f)^(3/2))).*((16*(B^2)/(SummationOfZz0(z,f,n)*F(z,f)))^((SummationOfZz0(z,f,n).*(2/B))^0.5))*exp((-1).*(2^(5/2)).*(B^(3/2))/(3*F(z,f))+(SummationOfZz0(z,f,n).*((2/B)^(0.5)))/3+(2^(5/2).*(term(z,f)^(3/2)))/(3*F(z,f)));
                        Ru(z,f)=(R(z,f)/u(z,f));
                    end
                end
                for f=1:1:Fsize
                     Pt(1,f)=1-exp((-1)*z0_resolution*sum(Ru(int32(SummationOfZc(20,f,n)+1):1:zsize,f)));
                    SummationPt(n,f)=Pt(1,f);
                end
            end
            for i=1:1:(size(SummationPt,1)-1)
%                 subplot(2,2,SubNumber);
                if i==1
                    format='-.k';
                elseif i==2
                    format='-k';
                elseif i==3
                    format=':k';
                elseif i==4
                    format='--k';
                elseif i==5
                    format='-ok';
                else
                    format='-xk';
                end
                
                semilogy(10:F_resolution:200,(SummationPt(i,10/F_resolution:Fsize)-SummationPt(i+1,10/F_resolution:Fsize)),format,'linewidth',3);
                xlabel('Field (V/nm)','fontsize',25);
                
                text(ConstantTable(e,4)*0.8,0.0001,ElementName,'fontsize',25,'BackgroundColor',[1 1 1]);
                axis([0 ConstantTable(e,4) 0.0001 1]);
                ax1=gca;
                set(ax1,'xaxislocation','top','linewidth',2.5,'fontsize',25);
                box off;
                set(gca,'ytick',[0.0001 0.001 0.01 0.1 1],'TickLength',[0.025 0.025]);
                if ConstantTable(e,4)<100
                    Interval=20;
                else
                    Interval=40;
                end
                set(gca,'xtick',[0 Interval*1 Interval*2 Interval*3 Interval*4 Interval*5 Interval*6 Interval*7 Interval*8 Interval*9 Interval*10 Interval*11 Interval*12 Interval*13 Interval*14 Interval*15 Interval*16 Interval*17 Interval*18 Interval*19 Interval*20],'TickLength',[0.025 0.025]);

                hold on;  
            end
            SubNumber=SubNumber+1;
            if (SubNumber>4)
                SubNumber=1;
            
                SaveFileName=['Set' num2str(e) '.tif'];
                saveas(gcf,SaveFileName,'tif');
                close all;
                figure('Position',[1,1,1440,900])
            end

end









