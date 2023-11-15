open F,$ARGV[0];

while(<F>){
	if ($_=~/^>.*/){@a=split " ",$_; print "$a[0]\n";}
	else {print;}
}
