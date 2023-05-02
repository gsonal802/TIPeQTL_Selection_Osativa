#!/usr/bin/perl

use strict ;
use warnings ;
use lib "/home/sg7367/perl/lib/perl5/";
use Statistics::Descriptive  ;
use List::Util qw[min max];

my $window_size = 25;
my $sliding = 5;
my $limit = 0;
my $window = 0;
my %header ;
my %pi ;
my %sites ;
my %count;
my %header2 ;
my %AFS ;       #$AFS{window}{pop}{scaffold}{pos}
my %data ;
my %fixed ;
my %nonzero ;

my %GstPrime_hash ;
my %seg_sites_hash ;
my %all_seg_sites_hash ;
my %theta_pi_hash ;
my %theta_W_hash ;
my %Taj_D_hash ;

my $sum_total = 0 ;

# Enter here the names of the populations to analyse (alphabetical order)
my @pops_to_analyze = ($ARGV[0], $ARGV[1]);
my $pop0 = $pops_to_analyze[0] ;
my $pop1 = $pops_to_analyze[1] ;
my $popNB0 = $ARGV[2];
my $popNB1 = $ARGV[3];
my $workdir = $ARGV[4];
my $outdir = $ARGV[5];
my $filename = $ARGV[6];
open TESTOUT, ">${outdir}/testout.FST.$pop0.$pop1.txt" ;

my %nuc ;
# Enter number of chromosomes for each pop
$nuc{$pop0}=$popNB0;
$nuc{$pop1}=$popNB1;

foreach my $pop (@pops_to_analyze){
        $sum_total += $nuc{$pop} ;
        print TESTOUT $pop, "\t" , $nuc{$pop}, "\n" ;
}

print TESTOUT $sum_total, "\n" ;
close TESTOUT ;

# Input file in format: Scaffold /t Position /t Derived alleles counts in pop0 /t Derived alleles counts in pop1 
open IN, "<$workdir/$filename" ;

while(<IN>){
        chomp $_ ;
        my @line = split(/\t/, $_) ;
        # Buffer array with values of every line
        if($_ =~ m/Chr/){
                        # Collect pop names from header line
                        foreach my $index (2..$#line){
                        $header{$index} = $line[$index] ;
                        }
                }
                else{   
                        my $sum = 0 ;
                        foreach my $index (2..$#line){
                        #print $line[$index], "\t***\t";
                        $sum += $line[$index] ;
                        }
                        # initialize snp count for each scaffold
                        if (! exists $count{$line[0]}){
                                $count{$line[0]}=0;
                        }
                        # Only increment count when SNP is variable
                        if($sum>0 && $sum<$sum_total){
                                $count{$line[0]} ++ ;
                        }
                        $limit = min(int($window_size/$sliding), int((($count{$line[0]}-1)/$sliding)+1));
                        #print $limit; 
                        if($count{$line[0]} <= $window_size){
                                $window = 0;
                        }
                        else{
                                $window = int((($count{$line[0]}-1) - $window_size)/$sliding)+1;
                        }
                        for(my $i = 0; $i<$limit; $i++) {
                                # Add all snp positions of the window to sites array
                                push @{$sites{$line[0]}{$window}}, $line[1] ;
                                #print $count{$line[0]}, "\t", $window, "\n";
                                # Add all allele counts to AFS array
                                foreach my $index (2..$#line){
                                $AFS{$line[0]}{$window}{$header{$index}}{$line[1]} = $line[$index] ;
                                }
                                $window++;
                        }
        }
    }

close IN ;


foreach my $scaff(keys %AFS){
        foreach my $win (sort{$a<=>$b} keys %{$AFS{$scaff}}){
                foreach my $pop (@pops_to_analyze){
                        $seg_sites_hash{$scaff}{$win}{$pop} = seg_sites(\%{$AFS{$scaff}{$win}{$pop}}, $nuc{$pop}) ;
                }
        }
}

foreach my $scaff(keys %AFS){
        foreach my $win (sort{$a<=>$b} keys %{$AFS{$scaff}}){
                $all_seg_sites_hash{$scaff}{$win} = all_seg_sites(\%{$AFS{$scaff}{$win}}, %nuc);

        # print TESTOUT $scaff, "\t" , $win, "\t" , $all_seg_sites_hash{$scaff}{$win}, "\n" ;
                if($all_seg_sites_hash{$scaff}{$win}>=$window_size){
                        $GstPrime_hash{$scaff}{$win} = GstPrime(\%{$AFS{$scaff}{$win}}, $pop0, $pop1, $nuc{$pop0}, $nuc{$pop1}) ;
                        foreach my $pop (@pops_to_analyze){
                                $theta_pi_hash{$scaff}{$win}{$pop} = theta_pi(\%{$AFS{$scaff}{$win}{$pop}}, $nuc{$pop}) ;
                                $theta_W_hash{$scaff}{$win}{$pop} = theta_W(\%{$AFS{$scaff}{$win}{$pop}}, $nuc{$pop}) ;
                                $Taj_D_hash{$scaff}{$win}{$pop} = Taj_D($theta_pi_hash{$scaff}{$win}{$pop}, $theta_W_hash{$scaff}{$win}{$pop}, $seg_sites_hash{$scaff}{$win}{$pop}, $nuc{$pop}) ;
                        }
                }
                else{
                        $GstPrime_hash{$scaff}{$win} = GstPrime(\%{$AFS{$scaff}{$win}}, $pop0, $pop1, $nuc{$pop0}, $nuc{$pop1}) ;
                        foreach my $pop (@pops_to_analyze){
                                $theta_pi_hash{$scaff}{$win}{$pop} = theta_pi(\%{$AFS{$scaff}{$win}{$pop}}, $nuc{$pop}) ;
                                if($seg_sites_hash{$scaff}{$win}{$pop} >= $window_size){
                                	$theta_W_hash{$scaff}{$win}{$pop} = theta_W(\%{$AFS{$scaff}{$win}{$pop}}, $nuc{$pop}) ;
                                	$Taj_D_hash{$scaff}{$win}{$pop} = Taj_D($theta_pi_hash{$scaff}{$win}{$pop}, $theta_W_hash{$scaff}{$win}{$pop}, $seg_sites_hash{$scaff}{$win}{$pop}, $nuc{$pop}) ;
                                }
                                else{
                     	       		$theta_W_hash{$scaff}{$win}{$pop} = 'NA';
                                	$Taj_D_hash{$scaff}{$win}{$pop} = 'NA';
                                }
                        }
                }
        }
}


my @tmp = split(/\./, $filename);
my $nameid = $tmp[0];
open OUT, ">${outdir}/Selstats.${nameid}.${window_size}SNPwin.txt" ;
print OUT "scaffold", "\t", "start", "\t", "end", "\t","seg_sites", "\t", "GSTprime", "\t", "pi_pr", "\t", "pi_ab", "\t", "TajD_pr","\t", "TajD_ab","\t";

foreach my $scaff(keys %AFS){
        foreach my $win (sort{$a<=>$b} keys %{$AFS{$scaff}}){
                print OUT "\n" ;
                print OUT $scaff, "\t" ;
                my $num_sites = 0 ;
                # $# calls the last index of the array
                print OUT ${$sites{$scaff}{$win}}[0], "\t", ${$sites{$scaff}{$win}}[$#{$sites{$scaff}{$win}}], "\t" ;
                $num_sites = ${$sites{$scaff}{$win}}[$#{$sites{$scaff}{$win}}]-${$sites{$scaff}{$win}}[0] ;
                print OUT $all_seg_sites_hash{$scaff}{$win}, "\t"  ;
                print OUT $GstPrime_hash{$scaff}{$win}, "\t"  ;
                foreach my $pop (@pops_to_analyze){
                print OUT sprintf("%.3e",$theta_pi_hash{$scaff}{$win}{$pop}/($num_sites)), "\t" ;
                }
                foreach my $pop (@pops_to_analyze){
                        print OUT $Taj_D_hash{$scaff}{$win}{$pop}, "\t" ;
                }
        }
}

close OUT;
exit ;

########### SUBROUTINES

# for calculations, see "Statistical Tests for Detecting Positive Selection by Utilizing High-Frequency Variants"


sub bc {
        my ($n,$k) = @_;
        my $r=1;
        $r*=$n/($n-$k),$n--while$n>$k;
        $r;
}

sub GstPrime{
        my ($AFS, $pop1_index, $pop2_index, $n_pop1, $n_pop2) = @_ ;
        my $Gst_total = 0 ;
        my $count = 0 ;
        my $num_pops = 2 ;
        foreach my $pos ( keys ( %{$$AFS{$pop1_index}} )){
                my $AF1 = $$AFS{$pop1_index}{$pos}/$n_pop1 ;
                my $AF2 = $$AFS{$pop2_index}{$pos}/$n_pop2 ;
                #print $AF1, $AF2;
                if( ($AF1+$AF2)>0 && ($AF1+$AF2)<($num_pops) ){
                        $count ++ ;
                        my $Gst = 0 ;
                        my $Gst_max = 0 ;
                        my $p_bar = (($AF1 + $AF2)/2);
                        my $Ht = 2*($p_bar)*(1-$p_bar) ;
                        my $Hs =  (2*$AF1*(1-$AF1) + 2*$AF2*(1-$AF2))/2;
                        if( $Ht == 0 ){
                                $Gst = 0 ;
                        }else{
                                $Gst = (($Ht - $Hs)/$Ht) ;
                                $Gst_max = ($num_pops-1)*(1-$Hs)/($num_pops-1+$Hs);
                                $Gst = $Gst/$Gst_max ;
                        }
                        $Gst_total += $Gst ;
                }
        }
        return(sprintf("%.3f",$Gst_total/$count)) ;
}

sub total_sites{
        my ($AFS) = @_ ;
        my $count = 0 ;
        foreach my $pos ( keys ( %$AFS )){
                $count ++ ;
        }
        return($count) ;
}


sub all_seg_sites{
        my ($AFS, %nuc) = @_ ;
        my %seg_sites ;
        foreach my $pop ( keys %$AFS ){
                my $n = $nuc{$pop};
                foreach my $pos ( keys %{$$AFS{$pop}}){
                        if(($$AFS{$pop}{$pos} > 0) && ($$AFS{$pop}{$pos} < $n)){
                                $seg_sites{$pos} = 1 ;
                        }
                }       
        }
        my $S = keys %seg_sites;
        return($S) ;
}
sub seg_sites{
         my ($AFS, $n) = @_ ;
        my $S = 0 ;
        my $S_singleton = 0 ;
        foreach my $pos ( keys ( %$AFS )){
                if(($$AFS{$pos} > 0) && ($$AFS{$pos} < $n)){
                        $S += 1 ;
                }
                if($$AFS{$pos} == 1){
                        $S_singleton ++ ;
                }
        }       
        return($S) ;
}

sub theta_pi{
        my ($AFS, $n) = @_ ;
        my $theta_pi = 0 ;
        foreach my $pos ( keys ( %$AFS )){
                $theta_pi += (($$AFS{$pos})*($n-$$AFS{$pos}))/bc($n,2) ;
        }
        return($theta_pi) ;
}

sub theta_W{
        my ($AFS, $n) = @_ ;
        my $theta_W = 0 ;
        my $a = 0 ;
        foreach ( 1 .. ($n - 1) ) { 
                $a = $a + 1/$_ ; 
        }
        foreach my $pos ( keys ( %$AFS )){
                if(($$AFS{$pos} > 0) && ($$AFS{$pos} < $n)){
                        $theta_W += 1/$a ;
                }
        }
        return($theta_W) ;
}


sub Taj_D{
        my($theta_pi, $theta_W, $S, $n) = @_ ;
        my $tajD = 0 ;
        my $a1 = 0 ; 
        my $a2 = 0 ;
        foreach ( 1 .. ($n - 1) ) { 
                $a1 = $a1 + 1/$_ ; 
                $a2 = $a2 + 1/($_)**2 ;
        }
        my $e1 = (1/$a1)*(($n+1)/(3*($n-1)) - (1/$a1)) ;
        my $e2 = (1/($a1**2 + $a2))*((2*(($n**2) + $n + 3))/(9*$n*($n-1)) - ($n+2)/($n*$a1) + ($a2)/($a1**2)) ;
        if ( (sqrt($S*$e1 + $S*($S-1)*$e2)) > 0 ) { 
                $tajD = ($theta_pi - $theta_W) / (sqrt($S*$e1 + $S*($S-1)*$e2)) ; 
                return (sprintf("%.3f",$tajD)) ;
        }else{
                return ("NA") ;
        }
}
