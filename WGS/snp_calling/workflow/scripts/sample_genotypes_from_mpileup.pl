#!/usr/bin/env perl
# Assumes VCFs are sorted and only one chromosome (contig)

use strict;
use warnings;
use Getopt::Long;
use v5.14;      # so srand returns the seed
use VCF; # https://metacpan.org/release/SNKWATT/VCF-1.003/view/lib/VCF.pm
# use Data::Dumper;

### GLOBALS
# VCF Column definitions
use constant {
    CHROM  => 0,
    POS    => 1,
    ID     => 2,
    REF    => 3,
    ALT    => 4,
    QUAL   => 5,
    FILTER => 6,
    INFO   => 7,
    FORMAT => 8,
    SAMPLE => 9,
};

# PL order (See description of GL and PL fields:
# https://samtools.github.io/hts-specs/VCFv4.2.pdf)
our @pl_order = (   "0/0", "0/1", "1/1", 
                    "0/2", "1/2", "2/2",
                    "0/3", "1/3", "2/3", "3/3",
                    "0/4", "1/4", "2/4", "3/4", "4/4" );

our $out_fh;

### GET OPTIONS
my $in_file = "/dev/stdin";
my $out_file = "/dev/stdout";

GetOptions ("in_file=s"   => \$in_file,
            "out_file=s"  => \$out_file)
or die("Error in command line arguments\n");

### SUBROUTINES

## Get null allele index from alt string
sub get_allele_index {
    my ($alt_str, $a) = @_;
    my @alt = split(/,/, $alt_str);

    my $out_i = -1;
    for(my $i=0; $i <= $#alt; $i++){
        if($alt[$i] eq $a){
            $out_i = $i;
            last;
        }
    }
    return $out_i;
}

## Convert PL to genotype probabilities (assuming sum to 1)
sub prob_from_PL {
    my ($pl_str) = @_;
    my @pl = split(/,/, $pl_str);
    my @gp;
    my $sum = 0.0;
    for(my $i=0; $i <= $#pl; $i++){
        $gp[$i] = 10 ** (- $pl[$i] / 10);
        $sum += $gp[$i];
    }

    for(my $i=0; $i <= $#gp; $i++){
        $gp[$i] = $gp[$i] / $sum;
    }
    return @gp;
}

## Sample genotype (returns index of @pl_order)
sub sample_GT {
    my ($gp_ref) = @_;
    my @gp = @{$gp_ref};

    my $prob_sum = 0;
    my $prob_targ = rand(1);
    for(my $i=0; $i <= $#gp; $i++){
        $prob_sum += $gp[$i];
        if($prob_sum > $prob_targ){
            return $i;
        }
    }
    return $#gp; # Return genotype index
}


## Is the null allele in the genotype(by GT index)?
sub null_allele_in_GT {
    my ($null_idx, $gt_idx) = @_;
    if($gt_idx == 0){
        return 0;
    } else {
        if($null_idx < 0){
            return 0;
        } else {
            my $a_idx = $null_idx + 1;
            if($a_idx == 1){
                return($gt_idx == 1 || $gt_idx ==  2 || $gt_idx ==  4 || $gt_idx ==  7 || $gt_idx == 11);
            } elsif($a_idx == 2){
                return($gt_idx == 3 || $gt_idx ==  4 || $gt_idx ==  5 || $gt_idx ==  8 || $gt_idx == 12);
            } elsif($a_idx == 3){
                return($gt_idx == 6 || $gt_idx ==  7 || $gt_idx ==  8 || $gt_idx ==  9 || $gt_idx == 13);
            } elsif($a_idx == 4){
                return($gt_idx ==10 || $gt_idx == 11 || $gt_idx == 12 || $gt_idx == 13 || $gt_idx == 14);
            }
        }
    }
    return 0;
}

## Remove the an allele from ALT string
sub remove_allele {
    my ($alt_str, $rm_idx) = @_;
    my @alt = split(/,/, $alt_str);
    splice(@alt, $rm_idx, 1);

    return join(",", @alt);
}

## Adjust genotypes if null allele was removed
sub adjust_gt_after_remove {
    my ($gt_str, $rm_idx) = @_;
    if($rm_idx < 0){
        return $gt_str;
    } else {
        my @gt = split(/\//, $gt_str);
        for(my $i = 0; $i <= $#gt; $i++){
            if($gt[$i] == ($rm_idx + 1)){
                die("Tried removing allele $rm_idx but GT is $gt_str\n");
            } elsif($gt[$i] > ($rm_idx + 1)){
                $gt[$i]--;
            }
        }
        
        return join("/", @gt);
    }

}

### MAIN
# Get seed
my $seed = srand();

# Load vcf
my $vcf = VCF->new(file=>$in_file);
$vcf->parse_header();

# Open output files
if($out_file =~ /\.gz$/){
    open $out_fh, "| bgzip -c > $out_file", or die $!;
} else {
    open $out_fh, ">", $out_file, or die $!;
}

# Remove PL header
$vcf->remove_header_line(key=>"FORMAT", ID=>"PL");

# Remove * header
$vcf->remove_header_line(key=>"ALT", ID=>'*');

# Add GT header
$vcf->add_header_line({key=>"FORMAT", ID=>"GT", Number=>1, Type=>"String", Description=>"Genotype"});

# Add seed
$vcf->add_header_line({key=>"sampleGenotypes", ID=>"seed", value=>$seed});


# Print header to output
print $out_fh $vcf->format_header();

# Read through VCF
my $line = $vcf->next_data_array();
while($#{$line} >= SAMPLE){    
    # Get index of PL and AD format fields
    my $pl_idx = $vcf->get_tag_index(${$line}[FORMAT],'PL',':');
    my $ad_idx = $vcf->get_tag_index(${$line}[FORMAT],'AD',':');

    # Skip if no PL field
    if($pl_idx == -1){
        $line = $vcf->next_data_array(); # Get next line
        next;
    }
    # Get array of PL strings
    my $pl_arr = $vcf->get_sample_field($line, $pl_idx);
    
    # Null allele index
    my $null_idx = get_allele_index(${$line}[ALT], '<*>');

    # Array of AD strings
    my $ad_arr = $vcf->get_sample_field($line, $ad_idx);

    # For each PL string
    my $do_print = 1;
    for(my $i = SAMPLE; $i <= $#{$line}; $i++){
        # Get GT probability array
        my @gp_arr = prob_from_PL(${$pl_arr}[$i - SAMPLE]);
        my $gt_idx = sample_GT(\@gp_arr);

        # Skip if sampled genotype contains the null field
        if(null_allele_in_GT($null_idx, $gt_idx)){
            $do_print = 0;
            last;
        } else {
            # Replace PL field with new GT field
            ${$line}[FORMAT] = $vcf->replace_field(${$line}[FORMAT], 'GT', $pl_idx, ':');
            ${$line}[$i]     = $vcf->replace_field(${$line}[$i], adjust_gt_after_remove($pl_order[$gt_idx], $null_idx), $pl_idx, ':');

            # Remove null allele from AD unless all four alleles present
            # print STDERR "NULL IDX = $null_idx\n";
            if($null_idx >= 0){
                my @ad_vals = split(/,/, ${$ad_arr}[$i - SAMPLE]);
                splice(@ad_vals, $null_idx + 1, 1);
                ${$line}[$i] = $vcf->replace_field(${$line}[$i], join(',', @ad_vals), $ad_idx, ':');
            }
        }
    }

    if($do_print){
        if($null_idx >= 0){
            # Remove null allele from ALT list
            ${$line}[ALT] = remove_allele(${$line}[ALT], $null_idx);
        }
        if(${$line}[ALT] eq ""){
            ${$line}[ALT] = ".";
        }
        
        print $out_fh $vcf->format_line($line);
    }

    $line = $vcf->next_data_array(); # Get next line
}

# Close files
$vcf->close();
close($out_fh);

