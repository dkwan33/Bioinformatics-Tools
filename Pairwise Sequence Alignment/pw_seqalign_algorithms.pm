#!/bin/perl
use strict;
use warnings;
use Exporter;
package PWSeqAlign;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
$VERSION     = '1.00';
@ISA         = qw(Exporter);
@EXPORT      = qw(); #these are exported by default
@EXPORT_OK   = qw(global_align_redux global_align local_align print_matrix); #these CAN be exported

#Author: David Kwan

#FUNCTION global_align_redux
#This function utilizes the needleman-wunsch algorithm to align two strings
#INPUT: STRING A (LENGTH N), STRING B (LENGTH M)
#OUTPUT: List of: Min Distance between String A and B,
#                 String C and String D which is the optimal alignment between A and B.
sub global_align_redux {
  my $a = shift; #string a
  my $b = shift; #string b
  my $n = length $a; #length of string a
  my $m = length $b; #length of string b
  my $d = []; #distance matrix
  my $t = []; #transition matrix

  #build first row and column of matrix d and t
  for (my $i = 0; $i <= $n; $i++) { #string a is the rows
    $d->[$i][0] = $i;
    $t->[$i][0] = 1; #the leftmost column of the matrix always points to the top
  }
  for (my $j = 0; $j <= $m; $j++) { #string b is the columns
    $d->[0][$j] = $j;
    $t->[0][$j] = 2; #the topmost row of the matrix always points to the left
  }

  #build the rest of matrix d
  for (my $i = 1; $i <= $n; $i++) { #rows
    for (my $j = 1; $j <= $m; $j++) { #columns

      #get the aligning characters and do a string equality comparison
      my $substr_a = substr($a, ($i - 1), 1);
      my $substr_b = substr($b, ($j - 1), 1);
      my $diag_score = ($substr_a eq $substr_b) ? 0 : 1; #ternary operator to evaluate diagonal

      #penalty scores are all set to 1
      my $top = $d->[$i - 1][$j] + 1; #vertical (previous row)
      my $left = $d->[$i][$j - 1] + 1; #horizontal (previous column)
      my $diag = $d->[$i - 1][$j - 1] + $diag_score; #diagonal (previous row and column)

      #find the minimum
      if ($diag <= $top && $diag <= $left) { #diag is the smallest of the three
        $d->[$i][$j] = $diag;
        $t->[$i][$j] = 0;
      } elsif ($top <= $diag && $top <= $left) { #top is the smallest of the three
        $d->[$i][$j] = $top;
        $t->[$i][$j] = 1;
      } else { #left is the smallest of the three
        $d->[$i][$j] = $left;
        $t->[$i][$j] = 2;
      }
    }
  }

  #traceback through the transition matrix to draw your alignment
  #in the transition matrix, 0 means diagonal, 1 means up, and 2 means left
  my $i = $n; #rows to step through
  my $j = $m; #columns to step through
  my $string1 = ""; #strings that we will build
  my $string2 = "";

  while ($i > 0 || $j > 0) {
    if ($t->[$i][$j] == 0) { #if the transition matrix says diagonal
      $string1 = substr($a, $i-1, 1) . $string1;
      $string2 = substr($b, $j-1, 1) . $string2;
      $i--;
      $j--;
    } elsif ($t->[$i][$j] == 1) { #if the transition matrix says top
      $string1 = substr($a, $i-1, 1) . $string1;
      $string2 = "_" . $string2;
      $i--;
    } else {
      $string1 = "_" . $string1; #if the transition matrix says left
      $string2 = substr($b, $j-1, 1) . $string2;
      $j--;
    }
  }

  return ($d->[$n][$m], $string1, $string2);
}

#-------------------------------------------------------------------------------
#FUNCTION global_align
#This function utilizes the needleman-wunsch algorithm to align two strings
#INPUT: STRING A (LENGTH N), STRING B (LENGTH M), GAP PENALTY, SCORING MATRIX
#OUTPUT: List of: Min Distance between String A and B,
#                 String C and String D which is the optimal alignment between A and B.
sub global_align {
  my $a = shift; #string a
  my $b = shift; #string b
  my $gap = shift; #gap penalty
  my $scoring = shift; #scoring matrix
  my $n = length $a; #length of string a
  my $m = length $b; #length of string b
  my $d = []; #distance matrix
  my $t = []; #transition matrix

  #build first row and column of matrix d and t
  for (my $i = 0; $i <= $n; $i++) { #string a is the rows
    $d->[$i][0] = $i;
    $t->[$i][0] = 1; #the leftmost column of the matrix always points to the top
  }
  for (my $j = 0; $j <= $m; $j++) { #string b is the columns
    $d->[0][$j] = $j;
    $t->[0][$j] = 2; #the topmost row of the matrix always points to the left
  }

  #build the rest of matrix d
  for (my $i = 1; $i <= $n; $i++) { #rows
    for (my $j = 1; $j <= $m; $j++) { #columns

      #get the aligning characters and reference the matrix
      my $substr_a = substr($a, ($i - 1), 1);
      my $substr_b = substr($b, ($j - 1), 1);
      my $diag_score = _diag_score($substr_a, $substr_b, $scoring);

      #calculate the scores that will be compared in minimization
      my $top = $d->[$i - 1][$j] + $gap; #vertical (previous row)
      my $left = $d->[$i][$j - 1] + $gap; #horizontal (previous column)
      my $diag = $d->[$i - 1][$j - 1] + $diag_score; #diagonal (previous row and column)

      #find the minimum
      if ($diag <= $top && $diag <= $left) { #diag is the smallest of the three
        $d->[$i][$j] = $diag;
        $t->[$i][$j] = 0;
      } elsif ($top <= $diag && $top <= $left) { #top is the smallest of the three
        $d->[$i][$j] = $top;
        $t->[$i][$j] = 1;
      } else { #left is the smallest of the three
        $d->[$i][$j] = $left;
        $t->[$i][$j] = 2;
      }
    }
  }

  #traceback through the transition matrix to draw your alignment
  #in the transition matrix, 0 means diagonal, 1 means up, and 2 means left
  my $i = $n; #rows to step through
  my $j = $m; #columns to step through
  my $string1 = ""; #strings that we will build
  my $string2 = "";

  while ($i > 0 || $j > 0) {
    if ($t->[$i][$j] == 0) { #if the transition matrix says diagonal
      $string1 = substr($a, $i-1, 1) . $string1;
      $string2 = substr($b, $j-1, 1) . $string2;
      $i--;
      $j--;
    } elsif ($t->[$i][$j] == 1) { #if the transition matrix says top
      $string1 = substr($a, $i-1, 1) . $string1;
      $string2 = "_" . $string2;
      $i--;
    } else {
      $string1 = "_" . $string1; #if the transition matrix says left
      $string2 = substr($b, $j-1, 1) . $string2;
      $j--;
    }
  }

  return ($d->[$n][$m], $string1, $string2);
}

#-------------------------------------------------------------------------------
#FUNCTION local_align
#This function utilizes the smith-waterman algorithm to align two strings
#INPUT: STRING A (LENGTH N), STRING B (LENGTH M), GAP PENALTY, SCORING MATRIX
#OUTPUT: List of: Highest score value in the distance matrix
#                 String C and String D which is the optimal alignment between A and B.
#                 Offsets to the start of the aligned portion.
sub local_align {
  my $a = shift; #string a
  my $b = shift; #string b
  my $gap = shift; #gap penalty
  my $scoring = shift; #scoring matrix
  my $n = length $a; #length of string a
  my $m = length $b; #length of string b
  my $d = []; #distance matrix
  my $t = []; #transition matrix

  #build first row and column of matrix d and t
  for (my $i = 0; $i <= $n; $i++) { #string a is the rows
    $d->[$i][0] = 0;
    $t->[$i][0] = 3;
  }
  for (my $j = 0; $j <= $m; $j++) { #string b is the columns
    $d->[0][$j] = 0;
    $t->[0][$j] = 3;
  }

  #build the rest of matrix d
  for (my $i = 1; $i <= $n; $i++) { #rows
    for (my $j = 1; $j <= $m; $j++) { #columns

      #get the aligning characters and reference the matrix
      my $substr_a = substr($a, ($i - 1), 1);
      my $substr_b = substr($b, ($j - 1), 1);
      my $diag_score = _diag_score($substr_a, $substr_b, $scoring);

      #calculate the scores that will be compared in maximization
      my $top = $d->[$i - 1][$j] + $gap; #vertical (previous row)
      my $left = $d->[$i][$j - 1] + $gap; #horizontal (previous column)
      my $diag = $d->[$i - 1][$j - 1] + $diag_score; #diagonal (previous row and column)

      #find the maximum
      if ($diag <= 0 && $top <= 0 && $left <= 0) { #if 0 is the largest of the four
        $d->[$i][$j] = 0;
        $t->[$i][$j] = 3;
      } elsif ($diag >= $top && $diag >= $left) { #diag is the largest of the four
        $d->[$i][$j] = $diag;
        $t->[$i][$j] = 0;
      } elsif ($top >= $diag && $top >= $left) { #top is the largest of the four
        $d->[$i][$j] = $top;
        $t->[$i][$j] = 1;
      } else { #left is the largest of the four
        $d->[$i][$j] = $left;
        $t->[$i][$j] = 2;
      }
    }
  }

  #find the maximum value in the distance matrix, and the indices of this value
  my @max = max_matrix($d, $n, $m);

  #traceback through the transition matrix to draw your alignment
  #in the transition matrix, 0 means diagonal, 1 means up, 2 means left, and 3 means stop
  my $i = $max[1]; #rows to step through
  my $j = $max[2]; #columns to step through
  my $string1 = ""; #strings that we will build
  my $string2 = "";

  while ($t->[$i][$j] != 3) {
    if ($t->[$i][$j] == 0) { #if the transition matrix says diagonal
      $string1 = substr($a, $i-1, 1) . $string1;
      $string2 = substr($b, $j-1, 1) . $string2;
      $i--;
      $j--;
    } elsif ($t->[$i][$j] == 1) { #if the transition matrix says top
      $string1 = substr($a, $i-1, 1) . $string1;
      $string2 = "_" . $string2;
      $i--;
    } else {
      $string1 = "_" . $string1; #if the transition matrix says left
      $string2 = substr($b, $j-1, 1) . $string2;
      $j--;
    }
  }

  return ($max[0], $string1, $string2, $i, $j);
}

#-------------------------------------------------------------------------------
#WORK FUNCTION _diag_score
#This function utilizes a given scoring matrix to score the difference between two nucleotides
#the scoring matrix is expected to be for ACGT in that order
#INPUT: character a, character b, reference to scoring matrix
#OUTPUT: number score
sub _diag_score {
  my $a = shift;
  my $b = shift;
  my $matrix = shift;
  my $row; #row to look at
  my $column; #column to look at

  if ($a eq "A") {
    $row = 0;
  } elsif ($a eq "C") {
    $row = 1;
  } elsif ($a eq "G") {
    $row = 2;
  } else {
    $row = 3;
  }

  if ($b eq "A") {
    $column = 0;
  } elsif ($b eq "C") {
    $column = 1;
  } elsif ($b eq "G") {
    $column = 2;
  } else {
    $column = 3;
  }

  return $matrix->[$row][$column];
}

#-------------------------------------------------------------------------------
#FUNCTION max_matrix
#This function finds the maximum value in a matrix
#INPUT: Reference to a matrix, dimensions of matrix n x m
#OUTPUT: number, and the indices in the matrix at which this number is found
sub max_matrix {
  my $matrix = shift;
  my $max = $matrix->[0][0];
  my $n = shift;
  my $m = shift;
  my $index1; #the indices where max is
  my $index2;

  for (my $i = 0; $i <= $n; $i++){
    for (my $j = 0; $j <= $m; $j++){
      if ($matrix->[$i][$j] > $max){
        $max = $matrix->[$i][$j];
        $index1 = $i;
        $index2 = $j;
      }
    }
  }
  return ($max, $index1, $index2);
}

#-------------------------------------------------------------------------------
#FUNCTION min
#This function finds the minimum of a list of numbers
#INPUT: List of numbers
#OUTPUT: number
sub min {
  my $min = shift;
  foreach (@_) {
    $min = $_ if $_ < $min;
  }
  return $min;
}

#-------------------------------------------------------------------------------
#FUNCTION print_matrix
#This function prints a matrix
#INPUT: Matrix, Dimension x of matrix, Dimension y of matrix
#OUTPUT: none
sub print_matrix {
  my $D = shift;
  my $n = shift;
  my $m = shift;

  for (my $i = 0; $i <= $n; $i++){
    for (my $j = 0; $j <= $m; $j++){
      print $D->[$i][$j] . " ";
    }
    print "\n";
  }
}

#-------------------------------------------------------------------------------
1;
