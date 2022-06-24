function [Name] = AdjustedNames_Plotting(testName)

if(strcmp('LumiraDX (Anterior Nasal Swab)',testName))
    Name='LumiraDX (Anterior nasal swab)';
elseif(strcmp('CareStart (Anterior Nasal Swab - FDA)',testName))    
    Name='CareStart (Anterior nasal swab)';
elseif(strcmp('BinaxNOW (FDA)',testName))  
    Name='BinaxNOW';
elseif(strcmp('Liaison (Anterior Nasal Swab)',testName))
    Name='Liaison (Anterior nasal swab)';
elseif(strcmp('Sofia (FDA)',testName))
    Name='Sofia';
elseif(strcmp('Liaison (Nasalpharyngeal Swab)',testName))
    Name='Liaison (Nasopharyngeal Swab)';
elseif(strcmp('LumiraDX (Nasopharyngeal Swab)',testName))
    Name='LumiraDX (Nasopharyngeal Swab)';
elseif(strcmp('CareStart (Nasopharyngeal Swab)',testName))
    Name='CareStart (Nasopharyngeal Swab)';
elseif(strcmp('SCoV-2',testName))
    Name='SCoV-2 Ag Detect';
elseif(strcmp('Status COVID+Flu',testName))
    Name='Status COVID-19/Flu';    
else
    Name=testName;
end

