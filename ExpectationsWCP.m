function bipartiteExpectationsWCP=ExpectationsWCP(active,mu,eta,etad,pd,ed,pz)
    px = 1-pz;

    rawExpectations(:,:)=coherentSourceChannel(active,mu,eta,etad,pd,ed,px);

    %convert 16-D pattern to 5-D POVM
    %Patterns:
    %column: [H,V,+,-,None]
    %row: Alice POVM
    MappingH=zeros(16,1);
    MappingH(9)=1; %1000
    MappingH(13)=0.5; %1100
    MappingV=zeros(16,1);
    MappingV(5)=1; %0100
    MappingV(13)=0.5; %1100
    
    MappingD=zeros(16,1);
    MappingD(3)=1; %0010
    MappingD(4)=0.5; %0011
    
    MappingA=zeros(16,1);
    MappingA(2)=1; %0001
    MappingA(4)=0.5; %0011
    
    MappingN=zeros(16,1);
    MappingN(1)=1; %0000
    MappingN(6)=1; %0101
    MappingN(7)=1; %0110
    MappingN(8)=1; %0111
    MappingN(10)=1; %1001
    MappingN(11)=1; %1010
    MappingN(12)=1; %1011
    MappingN(14)=1; %1101
    MappingN(15)=1; %1110
    MappingN(16)=1; %1111
    
    Mapping=[MappingH,MappingV,MappingD,MappingA,MappingN];
    bipartiteExpectationsWCP=rawExpectations(:,:)*Mapping;
end

function expectations=coherentSourceChannel(active,mu,eta,etad,pd,ed,px)
    expectations=zeros(4,16);
    pz=1-px;
    t=eta*etad; %total transmittance
    Poisson=@(mu,n) exp(-mu)*mu^n/factorial(n);

    theta=asin(sqrt(ed));
    PL=sin(pi/4-theta)^2;
    PU=cos(pi/4-theta)^2;
    
    mapping_passive=[[pz*(1-ed),pz*ed,px*PU,px*PL];[pz*ed,pz*(1-ed),px*PL,px*PU];[pz*PL,pz*PU,px*(1-ed),px*ed];[pz*PU,pz*PL,px*ed,px*(1-ed)]];
    mapping_active=[[1-ed,ed,PU,PL];[ed,1-ed,PL,PU];[PL,PU,1-ed,ed];[PU,PL,ed,1-ed]];
    
    for input=1:4
        %iterating over each input state
        
        for output=1:16
            %iterating over each pattern
            a=index1to4(output-1); %detector event, 4 elements corresponding to [H,V,D,A]
            
            Ppattern=1;
            if(active==0)
                %passive basis choice
                for k=1:4
                    %iterating over each detector
                    Pclick=1-Poisson(mu*t*mapping_passive(input,k),0); 
                    Pclick=1-(1-Pclick)*(1-pd); %effect of dark count
                    if(a(k)==1)
                        Ppattern=Ppattern*Pclick;
                    elseif(a(k)==0)
                        Ppattern=Ppattern*(1-Pclick);
                    end
                end
            else
                %active basis choice
                PpatternZ=1;
                prob_activeZ=[1,1,0,0];
                for k=1:4
                    %iterating over each detector
                    Pclick=(1-Poisson(mu*t*mapping_active(input,k),0));
                    Pclick=1-(1-Pclick)*(1-pd); %effect of dark count
                    Pclick=prob_activeZ(k)*Pclick; %effect of basis choice (active)
                    if(a(k)==1)
                        PpatternZ=PpatternZ*Pclick;
                    elseif(a(k)==0)
                        PpatternZ=PpatternZ*(1-Pclick);
                    end
                end
                PpatternX=1;
                prob_activeX=[0,0,1,1];
                for k=1:4
                    %iterating over each detector
                    Pclick=(1-Poisson(mu*t*mapping_active(input,k),0));
                    Pclick=1-(1-Pclick)*(1-pd); %effect of dark count
                    Pclick=prob_activeX(k)*Pclick; %effect of basis choice (active)
                    if(a(k)==1)
                        PpatternX=PpatternX*Pclick;
                    elseif(a(k)==0)
                        PpatternX=PpatternX*(1-Pclick);
                    end
                end
                Ppattern=px*PpatternX+pz*PpatternZ;
            end
            expectations(input,output)=Ppattern;
        end
    end
end

function a=index1to4(index)
    a(1) = floor(index/8);
    index = mod(index,8);
    a(2) = floor(index/4);
    index = mod(index,4);
    a(3) = floor(index/2);
    index = mod(index,2);
    a(4) = index;
end