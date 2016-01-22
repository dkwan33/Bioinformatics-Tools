#!/bin/perl
use warnings;
use strict;
package ClusteringTree;

#Author: David Kwan

#This module contains the implementation of tree structure for ClusteringAlg

# Creates a new clustering tree leaf with a single object, with the given
# label
sub new {
  my $label = $_[0];
  return bless {left => undef, right => undef, distance => undef, label => $label}, "ClusteringTree";
}

# Joins two clustering tree objects into one larger clustering tree, at the
# given distance
sub join {
  my $A = $_[0];
  my $B = $_[1];
  my $distance = $_[2];

  return bless {left => $A, right => $B, distance => $distance, label => undef}, "ClusteringTree";
}

# Recursively prints a clustering tree, in preorder traversal. Prints in the
# format (distance (left subtree) (right subtree)); prints (label) for
# leaves
sub print {
  my $tree = $_[0];

  print "(";

  if (defined $tree->{distance}) {
    print "$tree->{distance} ";
    $tree->{left}->print();
    print " ";
    $tree->{right}->print();
  } else {
    print "$tree->{label}";
  }
  print ")";
}

1;
