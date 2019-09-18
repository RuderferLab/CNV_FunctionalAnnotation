#script to annotate new CNV
# takes CNV file as input requiring the following space separated columns: id chromosome start_position end_position DEL/DUP

# run: perl annotate-cnv.pl [CNV_file] > [output_file]
# example first 5 lines, assumes there is a header
# id chr start end type
# CNV1 18 1316158 1321355 DEL
# CNV2 11 11964415 11966592 DEL
# CNV3 6 17500226 17506043 DUP
# CNV4 5 12741496 12746419 DEL

# output includes all genes that are affected by CNV, if a CNV does not affect any gene it will not be in the output


#file that has exon information for each gene
my $gene_file = "Ensemblv75/Ensemblv75_merged_exons";

# enhancer-gene target file
my $enh_gene_file = "Ensemblv75/PEC-enh-gene-links.bed";

# TAD information
my $tad_file = "chromatin_structure/DER-18_TAD_adultbrain.bed";

# read in gene/exon information
my %ginfo = ();
my @gene_map = ();
open(gene,$gene_file);
my $fline = <gene>;
while(<gene>){
    chomp($_);
    my @line = split(" ",$_);
    my @chr = split(/;/,$line[1]);
    $chr[0] =~ s/chr//g;
    my @starts = split(/;/,$line[2]);
    my @ends = split(/;/,$line[3]);
    $ginfo{$line[0]} = "$chr[0] $starts[0] $ends[-1]";
    push(@{$gene_map[$chr[0]]},$_);
}


print "id chr start end type len gene exon pro enh intad\n";
open(in, $ARGV[0]);
my $fline = <in>;
while(<in>){
    chomp($_);
    my $sv = $_;
    my @line = split(" ",$_);
    my $id = $line[0];
    my $svchr = $line[1];
    my $svstart = $line[2];
    my $svend = $line[3];
    my $svtype = $line[4];
    if( $svtype ne "DEL" && $svtype ne "DUP"){
	next;
    }
    my $svlen = $svend - $svstart;

    my %enhancer = ();
    my %all_genes = ();
    my @tads = ();

    open(tad, $tad_file);
    while(<tad>){
	chomp($_);
	my @t = split("\t",$_);
	$t[0] =~ s/chr//g;
	if( $t[0] == $svchr && $t[1] < $svend && $t[2] > $svstart ){
	    push(@tads, $_);
	}
    }

    open(enh, $enh_gene_file);
    while(<enh>){
	chomp($_);
	my @e = split("\t",$_);
	$e[0] =~ s/chr//g;
	my $e_prop = 0;
	if( $e[0] == $svchr && $e[1] < $svend && $e[2] > $svstart ){
	    my $min = $svend;
	    if( $e[2] < $svend){
		$min = $e[2];
	    }
	    my $max = $svstart;
	    if( $e[1] > $svstart ){
		$max = $e[1];
	    }
	    my $overlap = $min - $max;
	    $e_len = $e[2] - $e[1];
	    $e_prop = $overlap / $e_len;
	    $enhancer{$e[4]} += $e_prop;
	    $all_genes{$e[4]} = 1;
	}
    }

    my %gene = ();
    my %pro = ();
    for(my $i = 0; $i < scalar @{$gene_map[$svchr]}; $i++ ){
	my @g = split(" ",$gene_map[$svchr][$i]);
	my @starts = split(/;/,$g[2]);
	my @ends = split(/;/,$g[3]);
	my $exon_prop = 0;
	my $len = 0;
	my $pro_end = $starts[0] - 1;
	my $pro_start = $pro_end - 2000;
	if( $g[4] =~ /-/){
	    $pro_start = $ends[-1] + 1;
	    $pro_end = $pro_start + 2000;
	}
	if($svend > $pro_start && $svstart < $pro_end){
	    my $min = $svend;
            if( $pro_end < $svend){
                $min = $pro_end;
            }
            my $max = $svstart;
            if( $pro_start > $svstart ){
                $max = $pro_start;
            }
            my $overlap = $min - $max;
            $p_len = 2000;
            $p_prop = $overlap / $p_len;

            $pro{$g[0]} = $p_prop;
            $all_genes{$g[0]} = 1;
	}
	for( my $j = 0; $j < scalar @starts; $j++ ){
	    if($svend > $starts[$j] && $svstart < $ends[$j]){
		my $min = $svend;
		if( $ends[$j] < $svend){
		    $min = $ends[$j];
		}
		my $max = $svstart;
		if( $starts[$j] > $svstart ){
		    $max = $starts[$j];
		}
		my $overlap = $min - ($max-1);
		my $exonlength = $ends[$j] - $starts[$j];
		$len += $exonlength;

		$exon_prop += $overlap;
		$all_genes{$g[0]} = 1;
	    }
	}
	if( $exon_prop > 0){
	    my $prop = $exon_prop/$g[5];
	    $gene{$g[0]} = $prop;
	}
    }

    while ( my ($key, $value) = each(%all_genes) ) {
	my $e = $enhancer{$key};
	if( !exists $enhancer{$key}){
	    $e = 0
	}
	my $g = $gene{$key};
	if( !exists $gene{$key}){
	    $g = 0;
	}
	my $p = $pro{$key};
	if( !exists $pro{$key}){
	    $p = 0;
	}
	my $intad = 0;
	for( my $i = 0; $i < scalar @tads; $i++ ){
	    my @t = split("\t",$tads[$i]);
	    my @gi = split(" ",$ginfo{$key});

	    if($gi[2] > $t[1] && $gi[1] < $t[2]){
		$intad = 1;
	    }
	}
        print "$id $svchr $svstart $svend $svtype $svlen $key $g $p $e $intad\n";
    }
}
