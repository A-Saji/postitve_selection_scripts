
die "\n\tperl\tclean_MSA_AA.pl\t<Multiple-alignment-file (STRING)>\t<periodicity (peptide:1, CDS:3)>\n\n" if ($#ARGV!=1);


#UPLOAD THE AA SEQUENCE ALIGNMENT (FORMAT: FASTA)

open F,$ARGV[0] or die "Cant open file $ARGV[0]\n";
$i=-1;
while(<F>){
	if ($_=~/^>.*/) {$i++; push @idarr,$_;}
	else {
		chomp;
		push @{$twod[$i]},split "",$_;
	}
}
close (F);
#FIND LENGTH OF ALIGNMENT

$len=@{$twod[$i]};
#print "$twod[1][0]\n";


#FIND THE POSITIONS IN AA-ALIGNMENT WHERE FREQUENCY OF GAPS IS UPTO 25%
$intv=$ARGV[1];
for($j=0;$j<$len;$j+=$intv){
	$col="";
	for($k=0;$k<=$i;$k++){
		$col.=$twod[$k][$j];
	}
	$cnt=$col=~s/-/-/g;
	if ($cnt<=0.25*($i+1)){ 
		push @pos,$j;
	} 
}
#print "\n";


for($l=0;$l<=$i;$l++){
	print $idarr[$l];
	$seq=join "",@{$twod[$l]};
	for ($m=0;$m<$#pos;$m++){
		print substr ($seq,$pos[$m],$intv);
	}
	print "\n";
}



#open F,$ARGV[0] or die "Cant open file $ARGV[0]\n";
#while(<F>){
#        if ($_=~/^>.*/) {print }
#        else {
#                chomp;
#                for ($l=0; $l<=$#pos; $l++){
#			print substr ($_,$pos[$l],$intv);
#		}
#		print "\n";
#        }
#}
#close (F);



