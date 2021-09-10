    
function [cheng]=MCDZD(quik)
quik = strrep(quik,'U','');  % omit 'U'
    FF=quik;
    FFF=MCDtransform(FF);
    Fasta=MCDfen(FFF);
    Fasta1=MCD1D(Fasta);
    Fasta2=MCD2D(Fasta);
    Fasta3=MCD3D(Fasta);
    cheng=[Fasta1,Fasta2,Fasta3];
end