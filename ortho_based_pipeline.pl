
use warnings;
##Change the project path and app path according to your settings.
#make sure you have extract_seq_of_listedIDs.pl, translatorX perl scripts, clean_MSA_AA.pl, and Raxml installed in the project/scripts and app path.
#this code uses the orthogroups csv file as an input along with the gene sequences written here.
$path_proj="/home/angeosaji/scp/angeo_C4/Scpfolder";
$path_app="/home/angeosaji/scp/angeo_C4/Scpfolder";
#input the genes in the variable below
%c4gene2spc=('GRMZM2G039586','GRMZM2G039586','GRMZM2G471089','GRMZM2G471089')

open F,$ARGV[0] or die "Error: cant find Orthogroup file: $ARGV[0]\n";
while(<F>){
	foreach $id(keys%c4gene2spc){
		if ($_=~/$id/){
			chomp;
			my @cols=split " ",$_;
			my $genefam=$c4gene2spc{$id};
			open FO,">geneid_$genefam\.txt";
			for ($i=1;$i<=$#cols;$i++){
				$cols[$i]=~s/\.p$//;
				$cols[$i]=~s/\_P0[1-9]$//;
				#print "$cols[$i]\n";
				print FO "$cols[$i]\n";
			}
			close(FO);
			print "$genefam\n$_\n\n";
			system "perl $path_proj/scripts/extract_seq_of_listedIDs.pl geneid_$genefam.txt cds 300";
			system "perl $path_proj/scripts/extract_seq_of_listedIDs.pl geneid_$genefam.txt protein 100";
			system "perl $path_app/scripts/translatorx_vLocal.pl -i cds_300_$genefam.txt -o transAln_cds_300_$genefam";
			system "perl $path_proj/scripts/clean_MSA_AA.pl transAln_cds_300_$genefam.nt_ali.fasta 3 >transAln_cds_300_$genefam.nt_cleanali.fasta";
			system "perl $path_proj/scripts/clean_MSA_AA.pl transAln_cds_300_$genefam.aa_ali.fasta 1 >transAln_cds_300_$genefam.aa_cleanali.fasta";
			$r1=int(rand(100000));
			$r2=int(rand(100000));
			system "$path_app/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -f a -m PROTGAMMAAUTO -p $r1 -x $r2 -s transAln_cds_300_$genefam.aa_cleanali.fasta -# 200 -T 3 -n txt";
			open FO2,">tmp_codeml.ctl";
			print FO2 "seqfile = transAln_cds_300_$genefam.nt_cleanali.fasta\ntreefile = RAxML_bestTree.txt\noutfile = results_$genefam.txt\n\n";
			close (FO2);
			system "cat tmp_codeml.ctl $path_proj/scripts/codeml_fixedPart.txt >>codeml.ctl";
			system "$path_app/paml4.9j/bin/codeml";
			system "mkdir $genefam";
			system "mv *.* rst* rub* lnf $genefam";
		}
	}
}
