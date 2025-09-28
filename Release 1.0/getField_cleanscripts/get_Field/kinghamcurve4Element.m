function [KinghamCurves] = kinghamcurve4Element(ElementSymbol)
% Input Element is an Element symbol as  string, e.g. "Fe"




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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Find Element. Error if not found.
e = find(strcmp(ElementSymbol, ElementList));
    if isempty(e)
        error("I don't know this Element!");
    end


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
            n=n;
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
            
            % make a Table for the results
            colnames = ["Field", string(1:(size(SummationPt,1))) + "+"];           
            KinghamCurves = array2table(nan(length(10:F_resolution:200), length(colnames)), 'VariableNames', colnames);
            
            KinghamCurves.Field = (10:F_resolution:200)';
            
            %populate the table
            for k=1:1:(size(SummationPt,1)-1)
                KinghamCurves.(string(k) + "+") = (SummationPt(k,10/F_resolution:Fsize)-SummationPt(k+1,10/F_resolution:Fsize))';
            end
            %last column is just the highest chstate. Hope this is correct.
                KinghamCurves.(string(size(SummationPt,1))+"+") = (SummationPt(end,10/F_resolution:Fsize))';
            

            end









