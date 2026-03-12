#!/usr/bin/env perl
#
# Assumes columns:
# #[1]CHROM
# [2]POS
# [3]REF
# [4]ALT
# [5..n]<species>_<sample>:GT

use strict;
use warnings;
# use v5.10;

my ($out_format) = @ARGV;
$out_format //= "phylip"; # default 

$out_format = lc($out_format); # make lowercase

##### Functions #####
our %iupac_table = ( AA => "A",
                     TT => "T",
                     CC => "C",
                     GG => "G",
                     AG => "R",
                     GA => "R", 
                     CT => "Y",
                     TC => "Y",
                     GC => "S",
                     CG => "S",
                     AT => "W",
                     TA => "W",
                     GT => "K",
                     TG => "K",
                     AC => "M",
                     CA => "M");

# Get sequence character from a genotype
sub seqchar_from_gt{
    my ($ref, $alt, $gt) = @_;
    my @alleles = ($ref, split (/,/, $alt)); # Stuff the alleles into a single 

    if($gt eq "./." || $gt eq "." || $gt eq ""){
        return "N";
    } else {
        my $ntpair = $alleles[int(substr($gt, 0, 1))] . $alleles[int(substr($gt, 2, 1))]; # GT string to nt pair
        return $iupac_table{$ntpair};
    }
}

# Return true if ref and alt indicate that site is an indel (i.e., not snp)
sub is_snp{
    my ($ref, $alt) = @_;
    return !(length($ref) > 1 || (length($alt) > 1 && !($alt =~ /^[ATCG],[ATCG]$/ || $alt =~ /^[ATCG],[ATCG],[ATCG]$/)))
}


our $seqlen = 0; # Length of sequence of informative sites within window
our $windowsize = 0; # Size of current window
our %sequences;
our $samp_list = [];


# sub append_block{
#     my %cur_seqchar = %{$_[0]};

#     # Append characters to sequences in current block
#     foreach my $samp (@{$samp_list}){
#         $sequences{$samp} = $sequences{$samp} . $cur_seqchar{$samp}
#     }
#     $seqlen++; # Increase sequence length counter
# }

sub print_block{
    my ($out_format) = @_;

    if($seqlen > 0){
        if($out_format eq "phylip"){
            printf("%d %d\n", scalar(@{$samp_list}), $seqlen);
            foreach my $samp (@{$samp_list}){
                printf("%-45s%s\n", $samp, $sequences{$samp});
            }
            printf("\n");
        } elsif(($out_format eq "fa") || ($out_format eq "fasta")){
            foreach my $samp (@{$samp_list}){
                printf(">%s\n%s\n", $samp, $sequences{$samp});
            }
        }
        return 0;
    }
}

##### Process header #####
my @columns = split(' ', <STDIN>);

if($columns[0] ne "#[1]CHROM" || $columns[1] ne "[2]POS" || $columns[2] ne "[3]REF" || $columns[3] ne "[4]ALT"){
    # Exit if header is wrong 
    print STDERR "Column header must be in format:\n";
    print STDERR "[1]CHROM\n";
    print STDERR "[2]POS\n";
    print STDERR "[3]REF\n";
    print STDERR "[4]ALT\n";
    print STDERR "[5..n]<species>_<sample>:GT\n";
    exit 65;
}

# Initialize [sample -> column #] and [species -> sample] mappings
my %samp_to_col;
for(my $i = 4; $i <= $#columns; $i++){
    # Identify the columns and sample names associated with each column
    $columns[$i] =~ s/\[\d+\]//;
    $columns[$i] =~ s/:GT$//;
    
    $samp_to_col{$columns[$i]} = $i;
    push(@{$samp_list}, $columns[$i]);
}


##### Process rest of file #####
# my $num_sp_obs = 0; # Number of species observed
my $genotype = "";
# my %cur_seqchar; # Hash of current sequence characters

foreach my $samp (@{$samp_list}){
    $sequences{$samp} = "";
}

while(<STDIN>){
    my @columns = split;
    
    my $chr = $columns[0];
    my $pos = $columns[1];
    my $ref = $columns[2];
    my $alt = $columns[3];

    # Process the site
    # skip indels but not multi-allelics
    if(is_snp($ref, $alt)){
        foreach my $samp (@{$samp_list}){ # For each sample in species
            $genotype = $columns[$samp_to_col{$samp}];
            $sequences{$samp} = $sequences{$samp} . seqchar_from_gt($ref, $alt, $genotype); # Get the sequence character
        }
        $seqlen++; # Increase sequence length counter
    }
}

print_block($out_format);
