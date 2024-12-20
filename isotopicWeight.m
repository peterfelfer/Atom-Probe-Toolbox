function weight = isotopicWeight(Z,isotope)
% isotopicWeight gives you the weight of the element in the specific
% isotope state if input isotope is given. Otherwise all isotopes are listed
%
% weight = isotopicWeight(Z)
% weight = isotopicWeight(Z,isotope)
% 
% INPUT 
% Z:        atomic number
%
% isotope:  the isotope of the element (e.g., 56 for Fe), optional
%
% OUTPUT
% weight:   table with mass number, weight, and natural abundance of the
%           element's isotope(s)
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-N�rnberg

weight = []; % [isotope weight abundance]



switch Z
    case 1 %H
        weight = [1 1.01 99.985; 2 2.01 0.015];
        
    case 2 %He
        weight = [3 3.02 1.40E-04; 4 4.00 99.9999];
        
    case 3 %Li
        weight = [6 6.02 7.5; 7 7.02 92.5];

    case 4 %Be
        weight = [9 9.01 100];      
        
    case 5 %B
        weight = [10 10.01 19.9; 11 11.01 80.1];        
        
    case 6 %C
        weight = [12 12.00 98.89; 13 13.00 1.11];        
        
    case 7 %N
        weight = [14 14.00 99.634; 15 15.00 0.366];        
        
    case 8 %O
        weight = [16 15.995 99.762; 17 17.00 0.038; 18 18.00 0.2];        
        
    case 9 %F
        weight = [19 19.00 100];        
        
    case 10 %Ne
        weight = [20 19.99 90.48; 21 20.99 0.27; 22 21.99 9.25];        
        
    case 11 %Na
        weight = [23 22.99 100];       
        
    case 12 %Mg
        weight = [24 23.99 78.99; 25 24.99 10; 26 25.98 11.01];         
        
    case 13 %Al
        weight = [27 26.98 100];        
        
    case 14 %Si
        weight = [28 27.98 92.23; 29 28.98 4.67; 30 29.97 3.1];        
        
    case 15 %P
        weight = [31 30.97 100];        
        
    case 16 %S
        weight = [32 31.97 95.02; 33 32.97 0.75; 34 33.97 4.21; 36 35.97 0.02];        
        
    case 17 %Cl
        weight = [35 34.97 75.77; 37 36.97 24.23];        
        
    case 18 %Ar
        weight = [36 36.97 0.34; 38 37.96 0.06; 40 39.96 99.60];        
        
    case 19 %K
        weight = [39 38.96 93.258; 40 39.96 0.0117; 41 40.96 6.73];        
        
    case 20 %Ca
        weight = [40 39.96 96.941; 42 41.96 0.647; 43 42.96 0.135; 44 43.96 2.086; 46 45.95 0.004];       
        
    case 21 %Sc
        weight = [45 44.96 100];     
        
    case 22 %Ti
        weight = [46 45.95 8.25; 47 46.95 7.44; 48 47.95 73.72; 49 48.95 5.41; 50 49.95 5.18];       
        
    case 23 %V
        weight = [50 49.95 0.25; 51 50.94 99.75];
        
    case 24 %Cr
        weight = [50 49.95 4.345; 52 51.94 83.789; 53 52.94 9.501; 54 53.94 2.365];       
        
    case 25 %Mn
        weight = [55 54.94 100];       
        
    case 26 %Fe
        weight = [54 53.94 5.845; 56 55.94 91.754; 57 56.94 2.119; 58 57.93 0.282];       
        
    case 27 %Co
        weight = [59 58.93 100];       
        
    case 28 %Ni
        weight = [58 57.94 68.077; 60 59.93 26.223; 61 60.93 1.14; 62 61.93 3.634; 64 63.93 0.926];       
        
    case 29 %Cu
        weight = [63 62.93 69.17; 65 64.93 30.83];       
        
    case 30 %Zn
        weight = [64 63.93 48.6; 66 65.93 27.9; 67 66.93 4.1; 68 67.93 18.8; 70 69.93 0.6];       
        
    case 31 %Ga
        weight = [69 68.93 60.108; 71 70.93 39.892];       
        
    case 32 %Ge
        weight = [70 69.92 21.23; 72 71.92 27.66; 73 72.92 7.73; 74 73.92 35.94; 76 75.92 7.44];
        
    case 33 %As
        weight = [75 74.92 100];
        
    case 34 %Se
        weight = [74 73.92 0.89; 76 75.92 9.36; 77 76.92 7.63; 78 77.92 23.78; 80 79.92 49.61; 82 81.92 8.73];
        
    case 35 %Br
        weight = [79 78.92 50.69; 81 80.92 49.31];
        
    case 36 %Kr
        weight = [79 78.92 0.35; 80 79.92 2.25; 82 81.91 11.6; 83 82.91 11.5; 84 83.91 57; 86 85.91 17.3];        
        
    case 37 %Rb
        weight = [85 84.91 72.165; 87 86.91 27.835];
        
    case 38 %Sr
        weight = [84 83.91 0.56; 86 85.91 9.86; 87 86.91 7.02; 88 87.91 82.58];
        
    case 39 %Y
        weight = [89 88.91 100];
        
    case 40 %Zr
        weight = [90 89.91 51.45; 91 90.91 11.22; 92 91.91 17.15; 94 93.91 17.38; 96 95.91 2.8];
        
    case 41 %Nb
        weight = [93 92.91 100];
        
    case 42 %Mo
        weight = [92 91.91 14.84; 94 93.91 9.25; 95 94.91 15.92; 96 95.91 16.68; 97 96.91 9.55; 98 97.91 24.13; 100 99.91 9.63];
        
    case 44 %Ru
        weight = [96 95.91 5.52; 98 97.91 1.88; 99 98.91 12.7; 100 99.90 12.6; 101 100.91 17; 102 101.90 31.6; 104 103.91 18.7];
        
    case 45 %Rh
        weight = [103 102.91 100];        
        
    case 46 %Pd
        weight = [102 101.91 1.02; 104 103.90 11.14; 105 104.91 22.33; 106 105.90 27.33; 108 107.90 26.46; 110 109.91 11.72];
        
    case 47 %Ag
        weight = [107 106.91 51.839; 109 108.91 48.161];
        
    case 48 %Cd
        weight = [106 105.91 1.25; 108 107.90 0.89; 110 109.90 12.49; 111 110.90 12.8; 112 111.90 24.13; 113 112.90 12.22; 114 113.90 28.73; 116 115.90 7.49];
        
    case 49 %In
        weight = [113 112.90 4.29; 115 114.90 95.71]; 
        
    case 50 %Sn
        weight = [112 111.90 0.97; 114 113.90 0.65; 115 114.90 0.34; 116 115.90 14.54; 117 116.90 7.68; 118 117.90 24.22; 119 118.90 8.58; 120 119.90 32.59; 122 121.90 4.63; 124 123.91 5.79];
        
    case 51 %Sb
        weight = [121 120.90 57.21; 123 122.90 42.79];
        
    case 52 %Te
        weight = [120 119.90 0.096; 122 121.90 2.603; 123 122.90 0.908; 124 123.90 4.816; 125 124.90 7.139; 126 125.90 18.952; 128 127.90 31.687; 130 129.91 33.799];
        
    case 53 %I
        weight = [127 126.90 100];        
        
    case 54 %Xe
        weight = [124 123.91 0.1; 126 125.90 0.09; 128 127.90 1.91; 129 128.90 26.4; 130 129.90 4.1; 131 130.91 21.2; 132 131.90 26.9; 134 133.91 10.4; 136 135.91 8.9];
        
    case 55 %Cs
        weight = [133 132.91 100];
        
    case 56 %Ba
        weight = [130 129.91 0.106; 132 131.91 0.101; 134 133.90 2.417; 135 134.90 6.592; 136 135.90 7.854; 137 136.91 11.23; 138 137.91 71.7];
        
    case 58 %Ce
        weight = [136 135.91 0.19; 138 137.91 0.25; 140 139.91 88.48; 142 141.91 11.08];
        
    case 59 %Pr
        weight = [141 140.91 100];
        
    case 60 %Nd
        weight = [142 141.91 27.13; 143 142.91 12.18; 144 143.91 23.8; 145 144.91 8.3; 146 145.91 17.19; 148 147.92 5.76; 150 149.92 5.64];
        
    case 66 %Dy
        weight = [156 155.92 .056; 158 157.924 .095; 160 159.925 2.329; 161 160.9269 18.889; 162 161.9267 25.475; 163 162.928 24.896; 164 163.929 28.26];
    
    case 72 %Hf
        weight = [174 173.94 0.162; 176 175.94 5.206; 177 176.94 18.606; 178 177.94 27.297; 179 178.95 13.629; 180 179.95 35.1];
        
    case 73 %Ta
        weight = [180 179.95 0.012; 181 180.95 99.988];        
        
    case 74 %W
        weight = [180 179.95 0.12; 182 181.95 26.498; 183 182.95 14.314; 184 183.95 30.642; 186 185.95 28.426];
        
    case 75 %Re
        weight = [185 184.95 37.4; 187 186.96 62.6];
        
    case 76 %Os
        weight = [184 183.95 0.02; 186 185.95 1.58; 187 186.96 1.6; 188 187.96 13.3; 189 188.96 16.1; 190 189.96 26.4; 192 191.96 41];
        
    case 77 %Ir
        weight = [191 190.96 37.3; 193 192.96 63.7];
        
    case 78 %Pt
        weight = [190 189.96 0.01; 192 191.96 0.79; 194 193.96 32.9; 195 194.96 33.8; 196 195.96 25.3; 198 197.97 7.2];
        
    case 79 %Au
        weight = [197 196.97 100];
        
    case 80 %Hg
        weight = [196 195.97 0.15; 198 197.97 9.97; 199 198.97 16.87; 200 199.97 23.1; 201 200.97 13.18; 202 201.97 29.86; 204 203.97 6.87];
        
    case 82 %Pb
        weight = [204 203.97 1.4; 206 205.97 24.1; 207 206.98 22.1; 208 207.98 52.4];        
        
    case 83 %Bi
        weight = [209 208.98 100];
        
    case 90 %Th
        weight = [232 232.04 100];
        
    case 92 %U
        weight = [234 234.04 5.50E-03; 235 235.04 0.72; 238 238.05 99.27];
        
    otherwise
        
        disp(' isotope weight undefined!');
        weight = [0 0 0];
        
end




if exist('isotope','var')
    weight = weight(weight(:,1)==isotope,:);
    weight = array2table(weight);
    weight.Properties.VariableNames{1} = 'mass number';
    weight.Properties.VariableNames{2} = 'weight';
    weight.Properties.VariableNames{3} = 'natural abundance';
    weight.Properties.VariableUnits = {'', 'amu', '%'};
else

weight = array2table(weight);
weight.Properties.VariableNames{1} = 'mass number';
weight.Properties.VariableNames{2} = 'weight';
weight.Properties.VariableNames{3} = 'natural abundance';
weight.Properties.VariableUnits = {'', 'amu', '%'};
end

end
