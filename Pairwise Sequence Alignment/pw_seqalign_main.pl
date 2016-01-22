#!/usr/bin/perl
use strict;
use warnings;
use lib '.';
use pw_seqalign_algorithms qw( global_align_redux global_align local_align print_matrix );

#This program runs all the alignment functions in PWSeqAlign on some test cases.

my @result;

print "\n", "-"x40, "\n", "Testing Needleman-Wunsch without gap penalty or scoring matrix:\n";
my $testa = "GGGGA";
my $testb = "AAAAAAAAAATTGGGGATTT";
@result = PWSeqAlign::global_align_redux($testa, $testb);
print "Minimum Distance: ", shift @result, "\n";
print "Seq1: ", shift @result, "\n";
print "Seq2: ", shift @result, "\n";

print "\n", "-"x40, "\n", "Testing Needleman-Wunsch without gap penalty or scoring matrix:\n";
my $testa = "GGGGA";
my $testb = "GGGGA";
@result = PWSeqAlign::global_align_redux($testa, $testb);
print "Minimum Distance: ", shift @result, "\n";
print "Seq1: ", shift @result, "\n";
print "Seq2: ", shift @result, "\n";

print "\n", "-"x40, "\n", "Testing Needleman-Wunsch with gap penalty and scoring matrix:\n";
my $testa = "GGGGAC";
my $testb = "GGGGAT";
my $score = [[0, 1, 1, 0],[0, 1, 0, 1],[1, 2, 1, 1],[0, 1, 0, 0]];
@result = PWSeqAlign::global_align($testa, $testb, 5, $score);
print "Minimum Distance: ", shift @result, "\n";
print "Seq1: ", shift @result, "\n";
print "Seq2: ", shift @result, "\n";

print "\n", "-"x40, "\n", "Testing Needleman-Wunsch with gap penalty and scoring matrix:\n";
my $testa = "GGGGA";
my $testb = "AAAAAAAAAATTGGGGATTT";
my $score = [[0, 1, 1, 0],[0, 1, 0, 1],[1, 2, 1, 1],[0, 1, 0, 0]];
@result = PWSeqAlign::global_align($testa, $testb, 5, $score);
print "Minimum Distance: ", shift @result, "\n";
print "Seq1: ", shift @result, "\n";
print "Seq2: ", shift @result, "\n";


print "\n", "-"x40, "\n", "Testing Smith-Waterman with gap penalty and scoring matrix:\n";
my @result = PWSeqAlign::local_align($testa, $testb, -1, $score);
print "Max Alignment Score: ", shift @result, "\n";
print "Seq1: ", shift @result, "\n";
print "Seq2: ", shift @result, "\n";
print "Index1: ", shift @result, "\n";
print "Index2: ", shift @result, "\n";



# Blosum62 scoring matrix
print "\n", "-"x40, "\n", "Blosum62";
$" = "\n";
my $matrix2 = [[4,0,0,0],
             [0,9,-3,-1],
             [0,-3,6,-2],
             [0,-1,-2,5]];
print "\n", "-"x40, "\n", "SW-Alignment\n";
print "Max Alignment Score: ", "@{[&PWSeqAlign::local_align('ATTCGATCG', 'GTTGATGGC', 1, $matrix2)]}", "\n";
print "Max Alignment Score: ", "@{[&PWSeqAlign::local_align('AAA', 'AAAAAAGGGGG', 1, $matrix2)]}", "\n";
print "Max Alignment Score: ", "@{[&PWSeqAlign::local_align('ACACACTA', 'AGCACACA', 1, $matrix2)]}", "\n";

# Arbitrary matrix
# my $matrix3 = [[-4,0,0,0],
#              [0,-9,3,1],
#              [0,3,-6,2],
#              [0,1,2,-5]];

my $matrix3 = $matrix2;
print "\n", "-"x40, "\n", "NW-Alignment\n";
print "Minimum Distance: ", "@{[&PWSeqAlign::global_align('ATTCGATCG', 'ATTCGATCG', 1, $matrix3)]}", "\n";
print "Minimum Distance: ", "@{[&PWSeqAlign::global_align('GGGGA', 'GGGGA', 1, $matrix3)]}", "\n";
print "Minimum Distance: ", "@{[&PWSeqAlign::global_align('AAAAAAA', 'AAAAAAA', 1, $matrix3)]}", "\n";
print "Minimum Distance: ", "@{[&PWSeqAlign::global_align('ATTCGATCG', 'GTTGATGGC', 1, $matrix3)]}", "\n";
print "Minimum Distance: ", "@{[&PWSeqAlign::global_align('AAA', 'AAAAAAGGGGG', 1, $matrix3)]}", "\n";
print "Minimum Distance: ", "@{[&PWSeqAlign::global_align('ACACACTA', 'AGCACACA', 1, $matrix3)]}", "\n";
