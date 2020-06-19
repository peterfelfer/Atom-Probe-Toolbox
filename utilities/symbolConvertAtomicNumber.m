function symnum = symbolConvertAtomicNumber(input)
% symbolConvertAtomicNumber transform the atomic number of an element to his symbol and vice versa the
% symbol to the atomic number Z
% 
% symnum = symbolConvertAtomicNumber(input)
% 
% INPUT
% input:  number or char that is either the atomic number or the symbol of
%         the element
%
% OUTPUT
% symnum: symbol of the element or atomic number of the element

if exist('input','var')
    
    %% conversion from element symbol to atomic number
    if ischar(input)==1
        switch input
        case 'H'
            symnum = 1;

        case 'He'
            symnum = 2;

        case 'Li'
            symnum = 3;

        case 'Be'
            symnum = 4;

        case 'B'
            symnum = 5;

        case 'C' 
            symnum = 6;

        case 'N' 
            symnum = 7;

        case 'O'
            symnum = 8;

        case 'F'
            symnum = 9;

        case 'Ne'
            symnum = 10;

        case 'Na' 
            symnum = 11;

        case 'Mg'
            symnum = 12;

        case 'Al'
            symnum = 13;

        case 'Si' 
            symnum = 14;

        case 'P' 
            symnum = 15;

        case 'S'
            symnum = 16;

        case 'Cl'
            symnum = 17;

        case 'Ar'
            symnum = 18;

        case 'K' 
            symnum = 19;

        case 'Ca'
            symnum = 20;

        case 'Sc' 
            symnum = 21;

        case 'Ti' 
            symnum = 22;

        case 'V' 
            symnum = 23;

        case 'Cr' 
            symnum = 24;

        case 'Mn' 
            symnum = 25;

        case 'Fe' 
            symnum = 26;

        case 'Co' 
            symnum = 27;

        case 'Ni' 
            symnum = 28;

        case 'Cu' 
            symnum = 29;

        case 'Zn' 
            symnum = 30;

        case 'Ga' 
            symnum = 31;

        case 'Ge' 
            symnum = 32;

        case 'As' 
            symnum = 33;

        case 'Se'
            symnum = 34;

        case 'Br' 
            symnum = 35;

        case 'Kr' 
            symnum = 36;

        case 'Rb' 
            symnum = 37;

        case 'Sr' 
            symnum = 38;

        case 'Y' 
            symnum = 39;

        case 'Zr' 
            symnum = 40;

        case 'Nb' 
            symnum = 41;

        case 'Mo' 
            symnum = 42;

        case 'Ru'
            symnum = 44;

        case 'Rh' 
            symnum = 45;

        case 'Pd'
            symnum = 46;

        case 'Ag' 
            symnum = 47;

        case 'Cd'
            symnum = 48;

        case 'In'
            symnum = 49;

        case 'Sn'
            symnum = 50;

        case 'Sb' 
            symnum = 51;

        case 'Te'
            symnum = 52;

        case 'I'
            symnum = 53;

        case 'Xe'
            symnum = 54;

        case 'Cs'
            symnum = 55;

        case 'Ba'
            symnum = 56;

        case 'La'
            symnum = 57;

        case 'Ce'
            symnum = 58;

        case 'Pr'
            symnum = 59;

        case 'Nd'
            symnum = 60;

        case 'Pm'
            symnum = 61;

        case 'Sm'
            symnum = 62;

        case 'Eu'
            symnum = 63;

        case 'Gd'
            symnum = 64;

        case 'Tb'
            symnum = 65;

        case 'Dy'
            symnum = 66;

        case 'Ho'
            symnum = 67;

        case 'Er'
            symnum = 68;

        case 'Tm'
            symnum = 69;

        case 'Yb'
            symnum = 70;

        case 'Lu'
            symnum = 71;

        case 'Hf'
            symnum = 72;

        case 'Ta'
            symnum = 73;

        case 'W'
            symnum = 74;

        case 'Re'
            symnum = 75;

        case 'Os'
            symnum = 76;

        case 'Ir'
            symnum = 77;

        case 'Pt'
            symnum = 78;

        case 'Au'
            symnum = 79;

        case 'Hg'
            symnum = 80;

        case 'Tl'
            symnum = 81;

        case 'Pb'
            symnum = 82;

        case 'Bi'
            symnum = 83;

        case 'Po'
            symnum = 84;

        case 'At'
            symnum = 85;

        case 'Rn'
            symnum = 86;

        case 'Fr'
            symnum = 87;

        case 'Ra'
            symnum = 88;

        case 'Ac'
            symnum = 89;

        case 'Th'
            symnum = 90;

        case 'Pa'
            symnum = 91;

        case 'U'
            symnum = 92;

        case 'Np'
            symnum = 93;

        case 'Pu'
            symnum = 94;

        case 'Am'
            symnum = 95;

        case 'Cm'
            symnum = 96;

        case 'Bk'
            symnum = 97;

        case 'Cf'
            symnum = 98;

        case 'Es'
            symnum = 99;

        case 'Fm'
            symnum = 100;

        case 'Md'
            symnum = 101;

        case 'No'
            symnum = 102;

        case 'Lr'
            symnum = 103;

        case 'Rf'
            symnum = 104;

        case 'Db'
            symnum = 105;

        case 'Sg'
            symnum = 106;

        case 'Bh'
            symnum = 107;

        case 'Hs'
            symnum = 108;

        case 'Mt'
            symnum = 109;

        case 'Ds'
            symnum = 110;

        case 'Rg'
            symnum = 111;

        case 'Uub'
            symnum = 112;

        case 'Uzt'
            symnum = 113;

        case 'Uuq'
            symnum = 114;

        case 'Uup'
            symnum = 115;

        case 'Uuh'
            symnum = 116;

        case 'Uus'
            symnum = 117;

        case 'Uuo'
            symnum = 118;
            otherwise
        symnum = NaN;
        end
        
        %% conversion from atomic number to element symbol
    elseif ischar(input)==0
        switch input
        case 1 
            symnum = 'H';

        case 2 
            symnum = 'He';

        case 3
            symnum = 'Li';

        case 4
            symnum = 'Be';

        case 5 
            symnum = 'B';

        case 6 
            symnum = 'C';

        case 7 
            symnum = 'N';

        case 8 
            symnum = 'O';

        case 9 
            symnum = 'F';

        case 10 
            symnum = 'Ne';

        case 11 
            symnum = 'Na';

        case 12 
            symnum = 'Mg';

        case 13 
            symnum = 'Al';

        case 14 
            symnum = 'Si';

        case 15 
            symnum = 'P';

        case 16 
            symnum = 'S';

        case 17 
            symnum = 'Cl';

        case 18 
            symnum = 'Ar';

        case 19 
            symnum = 'K';

        case 20 
            symnum = 'Ca';

        case 21 
            symnum = 'Sc';

        case 22 
            symnum = 'Ti';

        case 23 
            symnum = 'V';

        case 24 
            symnum = 'Cr';

        case 25 
            symnum = 'Mn';

        case 26 
            symnum = 'Fe';

        case 27 
            symnum = 'Co';

        case 28 
            symnum = 'Ni';

        case 29 
            symnum = 'Cu';

        case 30 
            symnum = 'Zn';

        case 31 
            symnum = 'Ga';

        case 32 
            symnum = 'Ge';

        case 33 
            symnum = 'As';

        case 34 
            symnum = 'Se';

        case 35 
            symnum = 'Br';

        case 36 
            symnum = 'Kr';

        case 37 
            symnum = 'Rb';

        case 38 
            symnum = 'Sr';

        case 39 
            symnum = 'Y';

        case 40 
            symnum = 'Zr';

        case 41 
            symnum = 'Nb';

        case 42 
            symnum = 'Mo';

        case 43
            symnum = 'Tc';

        case 44 
            symnum = 'Ru';

        case 45 
            symnum = 'Rh';

        case 46 
            symnum = 'Pd';

        case 47 
            symnum = 'Ag';

        case 48 
            symnum = 'Cd';

        case 49 
            symnum = 'In';

        case 50 
            symnum = 'Sn';

        case 51 
            symnum = 'Sb';

        case 52 
            symnum = 'Te';

        case 53 
            symnum = 'I';

        case 54 
            symnum = 'Xe';

        case 55 
            symnum = 'Cs';        

        case 56 
            symnum = 'Ba';

        case 57 
            symnum = 'La';

        case 58 
            symnum = 'Ce';

        case 59 
            symnum = 'Pr';

        case 60 
            symnum = 'Nd';

        case 61
            symnum = 'Pm';

        case 62
            symnum = 'Sm';

        case 63
            symnum = 'Eu';

        case 64
            symnum = 'Gd';

        case 65
            symnum = 'Tb';

        case 66
            symnum = 'Dy';

        case 67
            symnum = 'Ho';

        case 68
            symnum = 'Er';

        case 69
            symnum = 'Tm';

        case 70
            symnum = 'Yb';

        case 71
            symnum = 'Lu';

        case 72 
            symnum = 'Hf';

        case 73 
            symnum = 'Ta';

        case 74 
            symnum = 'W';

        case 75 
            symnum = 'Re';

        case 76 
            symnum = 'Os';

        case 77 
            symnum = 'Ir';

        case 78 
            symnum = 'Pt';

        case 79 
            symnum = 'Au';

        case 80 
            symnum = 'Hg';

        case 81
            symnum = 'Tl';

        case 82 
            symnum = 'Pb';

        case 83 
            symnum = 'Bi';

        case 84
            symnum = 'Po';

        case 85
            symnum = 'At';

        case 86
            symnum = 'Rn';

        case 87
            symnum = 'Fr';

        case 88
            symnum = 'Ra';

        case 89
            symnum = 'Ac';

        case 90 
            symnum = 'Th';

        case 91
            symnum = 'Pa';

        case 92 
            symnum = 'U';

        case 93
            symnum = 'Np';

        case 94
            symnum = 'Pu';

        case 95
            symnum = 'Am';

        case 96
            symnum = 'Cm';

        case 97
            symnum = 'Bk';

        case 98
            symnum = 'Cf';

        case 99
            symnum = 'Es';

        case 100
            symnum = 'Fm';

        case 101
            symnum = 'Md';

        case 102
            symnum = 'No';

        case 103
            symnum = 'Lr';

        case 104
            symnum = 'Rf';

        case 105
            symnum = 'Db';

        case 106
            symnum = 'Sg';

        case 107
            symnum = 'Bh';

        case 108
            symnum = 'Hs';

        case 109
            symnum = 'Mt';

        case 110
            symnum = 'Ds';

        case 111
            symnum = 'Rg';

        case 112
            symnum = 'Uub';

        case 113
            symnum = 'Uut';

        case 114
            symnum = 'Uuq';

        case 115
            symnum = 'Uup';

        case 116
            symnum = 'Uuh';

        case 117
            symnum = 'Uus';

        case 118
            symnum = 'Uuo'; 
        otherwise
        symnum = 'undefined';
        end
    end
else
    disp(' no input!');
end
