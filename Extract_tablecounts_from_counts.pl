# Cecile Pereira
# 10 Oct 2016
# Extract table from htseq count results
# perl Extract_tablecounts_from_counts.pl ../HTSEQC_150_Trim/ > ../table_reads_counts.csv

$DATA=$ARGV[0];#directory with htseqcount result
chdir($DATA);

$sample='';
%genesample=();
%samplestot=();
@counts=glob("*counts.txt");
foreach $f(@counts){
	open(IN,$f);
	if($f=~/ah_(.+)_counts|.txt/){
		$sample=$1;
		$samplestot{$sample}='';
	}
	while(<IN>){
			chomp;
			@tmp=split(/\s/,$_);
			$genesample{$tmp[0]}{$sample}=$tmp[1];
	}	
	close(IN);
}
@tablesample=keys(%samplestot);
$dernsample=$tablesample[$#tablesample];
$namesamples=join("\t",keys(%samplestot));
print "$namesamples\n";
foreach $gene(keys(%genesample)){
  print "$gene\t";
  foreach $i (keys(%samplestot)){
    if(exists($genesample{$gene}{$i})){
			if($genesample{$gene}{$i} eq "0"){
				print "0.00";
			}
			else{
      		print "$genesample{$gene}{$i}";
			}
    }
    else{
      print "0.00";
    }
    unless($i eq $dernsample){
    	#print "\n $i $dernsample dans le unless\n";
	    print "\t";
    }	
  }
  print "\n";
}

