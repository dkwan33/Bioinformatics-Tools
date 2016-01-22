#!/bin/perl
use strict;
use warnings;
use lib '.';
use clust_tree;
package ClusteringAlg;


#Author: David Kwan

#This module contains the implementation of four hierarchiccal clustering algorithms.
#UPGMA, WPGMA, Single Linkage, and Complete Linkage.

#-------------------------------------------------------------------------------
#FUNCTION upgma
#This function performs the UPGMA clustering algorithm
#INPUT: D[n][n] (square matrix of distances between n objects)
#OUTPUT: clustering tree
sub upgma {
  my $D = shift;
  my @trees; #declare the array of clustering trees
  my @weight; #declare the array of weights
  my $chain = []; #declare array reference of stack of objects
  my $n = scalar @$D;

  #create new leaves and set their weights
  for (my $i = 0; $i < $n; $i++){
    $trees[$i] = ClusteringTree::new($i + 1);
    $weight[$i] = 1;
  }

  for (my $x = 0; $x < $n - 1; $x++) { #while number of objects > 1
    if (@$chain == 0) { #arbitrarily choose the first object
      push @$chain, 0;
    }

    #find the two closest objects
    my ($a, $b) = find_min_loop($D, $chain);

    if ($a > $b) {
      ($a, $b) = ($b, $a); #swap a and b, so that we keep the lower number as a
    }

    $trees[$a] = ClusteringTree::join($trees[$a], $trees[$b], $D->[$a][$b]);

    for (my $i = 0; $i < $n; $i++){
      next if($i == $a || $i == $b); #skip a and b
      #find the weighted average
      $D->[$a][$i] = ($D->[$a][$i] * $weight[$a] + $D->[$b][$i] * $weight[$b]) / ($weight[$a] + $weight[$b]);
      $D->[$i][$a] = $D->[$a][$i]; #keep the matrix symmetric
      $D->[$b][$i] = 9999999999; #set to infinity
      $D->[$i][$b] = $D->[$b][$i]; #set to infinnity
    }

    $D->[$a][$b] = 9999999999;
    $D->[$b][$a] = $D->[$a][$b];

    $weight[$a] += $weight[$b];
  }

  return $trees[0]; #the full clustering tree has been collapsed into the first element
}

#-------------------------------------------------------------------------------
#FUNCTION wpgma
#This function performs the WPGMA clustering algorithm
#INPUT: D[n][n] (square matrix of distances between n objects)
#OUTPUT: clustering tree
sub wpgma {
  my $D = shift;
  my @trees; #declare the array of clustering trees
  my @weight; #declare the array of weights
  my $chain = []; #declare array reference of stack of objects
  my $n = scalar @$D;

  #create new leaves and set their weights
  for (my $i = 0; $i < $n; $i++){
    $trees[$i] = ClusteringTree::new($i + 1);
    $weight[$i] = 1;
  }

  for (my $x = 0; $x < $n - 1; $x++) { #while number of objects > 1
    if (@$chain == 0) {#arbitrarily choose the first object
      push @$chain, 0;
    }

    #find the two closest objects
    my ($a, $b) = find_min_loop($D, $chain);

    if ($a > $b) {
      ($a, $b) = ($b, $a); #swap a and b, so that we keep the lower number as a
    }

    $trees[$a] = ClusteringTree::join($trees[$a], $trees[$b], $D->[$a][$b]);

    for (my $i = 0; $i < $n; $i++){
      next if($i == $a || $i == $b); #skip a and b
      #find the average (not weighted)
      $D->[$a][$i] = ($D->[$a][$i] + $D->[$b][$i]) / 2;
      $D->[$i][$a] = $D->[$a][$i]; #keep the matrix symmetric
      $D->[$b][$i] = 9999999999; #set to infinity
      $D->[$i][$b] = $D->[$b][$i]; #set to infinnity
    }

    $D->[$a][$b] = 9999999999;
    $D->[$b][$a] = $D->[$a][$b];

    $weight[$a] += $weight[$b];
  }

  return $trees[0]; #the full clustering tree has been collapsed into the first element
}

#-------------------------------------------------------------------------------
#FUNCTION single_linkage
#This function performs the single linkage clustering algorithm
#INPUT: D[n][n] (square matrix of distances between n objects)
#OUTPUT: clustering tree
sub single_linkage {
  my $D = shift;
  my @trees; #declare the array of clustering trees
  my @weight; #declare the array of weights
  my $chain = []; #declare array reference of stack of objects
  my $n = scalar @$D;

  #create new leaves and set their weights
  for (my $i = 0; $i < $n; $i++){
    $trees[$i] = ClusteringTree::new($i + 1);
    $weight[$i] = 1;
  }

  for (my $x = 0; $x < $n - 1; $x++) { #while number of objects > 1
    if (@$chain == 0) { #arbitrarily choose the first object
      push @$chain, 0;
    }

    #find the two closest objects
    my ($a, $b) = find_min_loop($D, $chain);

    if ($a > $b) {
      ($a, $b) = ($b, $a); #swap a and b, so that we keep the lower number as a
    }

    $trees[$a] = ClusteringTree::join($trees[$a], $trees[$b], $D->[$a][$b]);

    for (my $i = 0; $i < $n; $i++){
      next if($i == $a || $i == $b); #skip a and b
      #find the minimum and use the minimum
      if ($D->[$a][$i] < $D->[$b][$i]) {
        $D->[$a][$i] = $D->[$a][$i];
      } else {
        $D->[$a][$i] = $D->[$b][$i];
      }
      $D->[$i][$a] = $D->[$a][$i]; #keep the matrix symmetric
      $D->[$b][$i] = 9999999999; #set to infinity
      $D->[$i][$b] = $D->[$b][$i]; #set to infinnity
    }

    $D->[$a][$b] = 9999999999;
    $D->[$b][$a] = $D->[$a][$b];

    $weight[$a] += $weight[$b];
  }

  return $trees[0]; #the full clustering tree has been collapsed into the first element
}

#-------------------------------------------------------------------------------
#FUNCTION complete_linkage
#This function performs the complete linkage clustering algorithm
#INPUT: D[n][n] (square matrix of distances between n objects)
#OUTPUT: clustering tree
sub complete_linkage {
  my $D = shift;
  my @trees; #declare the array of clustering trees
  my @weight; #declare the array of weights
  my $chain = []; #declare array reference of stack of objects
  my $n = scalar @$D;

  #create new leaves and set their weights
  for (my $i = 0; $i < $n; $i++){
    $trees[$i] = ClusteringTree::new($i + 1);
    $weight[$i] = 1;
  }

  for (my $x = 0; $x < $n - 1; $x++) { #while number of objects > 1
    if (@$chain == 0) { #arbitrarily choose the first object
      push @$chain, 0;
    }

    #find the two closest objects
    my ($a, $b) = find_min_loop($D, $chain);

    if ($a > $b) {
      ($a, $b) = ($b, $a); #swap a and b, so that we keep the lower number as a
    }

    $trees[$a] = ClusteringTree::join($trees[$a], $trees[$b], $D->[$a][$b]);

    for (my $i = 0; $i < $n; $i++){
      next if($i == $a || $i == $b); #skip a and b
      #find the maximum and use the maximum
      if ($D->[$a][$i] > $D->[$b][$i]) {
        $D->[$a][$i] = $D->[$a][$i];
      } else {
        $D->[$a][$i] = $D->[$b][$i];
      }
      $D->[$i][$a] = $D->[$a][$i]; #keep the matrix symmetric
      $D->[$b][$i] = 9999999999; #set to infinity
      $D->[$i][$b] = $D->[$b][$i]; #set to infinnity
    }

    $D->[$a][$b] = 9999999999;
    $D->[$b][$a] = $D->[$a][$b];

    $weight[$a] += $weight[$b];
  }

  return $trees[0]; #the full clustering tree has been collapsed into the first element
}

#-------------------------------------------------------------------------------
#FUNCTION: find_min_loop
#This function finds the two closest objects in a stack of objects
#INPUT: D[n][n] (distance matrix), chain (stack of objects)
#OUTPUT: (a,b) (pair of objects such that the closest object to a is b and vice versa)
sub find_min_loop {
  my $D = shift;
  my $chain = shift;

  my $cur_obj = $chain->[-1]; #get the current object at the end of the chain

  while (1) {
    my $min_dist = 99999999999; #set to infinity
    my $near_obj; #this is the nearest object to cur_obj
    for (my $j = 0; $j < @$D; $j++){
      next if $j == $cur_obj; #skip cur_obj

      if ($D->[$cur_obj][$j] < $min_dist) { #if the distance is less than the current minimum distance
        $near_obj = $j; #this is the new closest object
        $min_dist = $D->[$cur_obj][$j]; #use the new minimum distance
      }
    }
    foreach my $element (@$chain) {
      if ($near_obj == $element) { #a loop is found
        pop @$chain; #merge the top object and its nearest neighbour
        pop @$chain;

        return ($cur_obj, $near_obj); #exit the find_min_loop function
      }
    }
    #no loop has been found, continue trying
    push @$chain, $near_obj;
    $cur_obj = $near_obj;
  }
}

1;
