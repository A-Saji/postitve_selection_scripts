use warnings;
#use strict;

die "\n\tUsage: perl extract_seq_of_listedIDs.pl\t<file_with_ID-List_in_1st_column (STR)>\t<molecule_type (cds/protein)>\t<minimum_length (INT)>\n\n" if ($#ARGV!=2);

#VARIABLES TO BE USED LATER	
	my %pref2spc=("LOC_Os","rice","Sobic","sorghum","GRMZM","maize","AT","Athaliana","Seita","Sitalica","Pavir","Pvirgatum","Sevir","Sviridis","Pahal","Phallii","Brast","Bstacei","Zm00008a","ZmaysPH207","Bradi","Bdistachyon","Chala","Chasmanthum","OEL","Dicanthelium");
	my %spc2file=("rice","Osativa_323_v7.0","sorghum","Sbicolor_313_v3.1","maize","Zmays_284_Ensembl-18_2010-01-MaizeSequence","Athaliana","Athaliana_167_TAIR10","Sitalica","Sitalica_312_v2.2","Pvirgatum","Pvirgatum_450_v4.1","Sviridis","Sviridis_500_v2.1","Phallii","Phallii_495_v3.1","Bstacei","Bstacei_316_v1.1","ZmaysPH207","ZmaysPH207_443_v1.1","Bdistachyon","Bdistachyon_314_v3.1","Dicanthelium","Dicanthelium1","Chasmanthum","Claxum_690_v1.1");
	my %spc2id=();
	#my $ref_file_path="/home/vivek/res_work/ref_genomes";
	my $ref_file_path="/home/angeosaji/C4_genes/SCP/ref_genomes";
	my $swi=0;
	my %id2seq=();
	my $outfile=$ARGV[0]; $outfile=~s/geneid_//;
#	open FOUT,">prot_$outfile" or die "\n\tError: Cant open output file: prot_$outfile\n\n";
#	open FOUT,">cds_$outfile" or die "\n\tError: Cant open output file: cds_$outfile\n\n";
	open FOUT,">$ARGV[1]\_$ARGV[2]\_$outfile" or die "\n\tError: Cant open output file: $ARGV[1]\_$ARGV[2]\_$outfile\n\n";
	my $len=0;
	my $seq=();

#OPEN AND UPLOAD IDs
	open FH,$ARGV[0] or die "Error: Cant locate file $ARGV[0]\n";
	while(<FH>){
		if ($_!~/^#/){
			chomp; 
			my @cols=split "\t",$_;
			$cols[0]=~s/,$//;
			foreach my $pref(keys%pref2spc){
				if ($cols[0]=~/$pref/){
					my $spc=$pref2spc{$pref};
					if (exists $spc2id{$spc}){push @{$spc2id{$spc}},$cols[0]; }
					else {$spc2id{$spc}[0]=$cols[0];}
				}
			} 
			#print "$cols[0]\t$pref\t$pref2spc{$pref}\n";
		}
	}
	close (FH);

#EXTRACT CDS SEQUENCES
	foreach my $spc(keys%spc2id){
		my $file=$spc2file{$spc};
		#$file.='.cds.fa';
#		$file.='.protein_primaryTranscriptOnly.fa';
#		$file.='.cds_primaryTranscriptOnly.fa';
		$file.="\.$ARGV[1]\_primaryTranscriptOnly.fa";
		print "$file\n";
		open FH,"$ref_file_path/$spc/$file" or die "Error: Cant open sequencefile: $ref_file_path/$spc/$file\n";

		while(<FH>){
			if (($swi==1)&&($_!~/^>.*/)) {$seq.=$_; chomp; $len+=length($_);}	#print FOUT $_;
			elsif ($_=~/^>.*/) {
				if ($swi==1) {print "$len\n"; if ($len>=$ARGV[2]) {print FOUT $seq;}}
#				if ($swi==1) {print "$len\n"; if ($len>=2500) {print FOUT $seq;}}
				$swi=0;
				#splice($_,0,1);
				my @cols=split /\s+/,$_;
				foreach my $id(@{$spc2id{$spc}}){
					foreach $col(@cols){
						if ($cols[0]=~/$id/) {print "$id\t$cols[0]\t"; $swi=1; $seq="\>$id\n"; $len=0; last;} #print FOUT $_;
						elsif ($col eq "polypeptide=$id") {print "$id\t$cols[0]\t"; $swi=1; $seq="\>$id\n"; $len=0; last;} #print FOUT $_;
						elsif ($col eq "locus=$id") {print "$id\t$cols[0]\t"; $swi=1; $seq="\>$id\n"; $len=0; last;} #print FOUT $_;
						#if ($cols[0]=~/$id/) {print "$id\t$cols[0]\t"; $swi=1; $seq="\>$id\n"; $len=0;} #print FOUT $_;
					}
				}
			}
		}

		close (FH);
		#print "$id2seq{'AT2G42600'}";
		#last;
	}
