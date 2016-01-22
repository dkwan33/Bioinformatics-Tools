#!/bin/perl
use strict;
use warnings;
use lib '.';
use clust_algorithms;
use clust_tree;

#This program runs all the clustering algorithms on some test cases.

#First few matrices here:http://www.southampton.ac.uk/~re1u06/teaching/upgma/
my $d1 = [[10,.19,.27,.08,.33,.18,.13],
	[.19,10,.31,.18,.36,.01,.13],
	[.27,.31,10,.26,.41,.32,.29],
	[.08,.18,.26,10,.31,.17,.14],
	[.33,.36,.41,.31,10,.35,.28],
	[.18,.01,.32,.17,.35,10,.12],
	[.13,.13,.29,.14,.28,.12,10]];

my $d2 = [
	[10,.19,.27,.08,.33,.18,.13],
	[.19,10,.31,.18,.36,.01,.13],
	[.27,.31,10,.26,.41,.32,.29],
	[.08,.18,.26,10,.31,.17,.14],
	[.33,.36,.41,.31,10,.35,.28],
	[.18,.01,.32,.17,.35,10,.12],
	[.13,.13,.29,.14,.28,.12,10]];

my $d3 = [
	[10,.19,.27,.08,.33,.18,.13],
	[.19,10,.31,.18,.36,.01,.13],
	[.27,.31,10,.26,.41,.32,.29],
	[.08,.18,.26,10,.31,.17,.14],
	[.33,.36,.41,.31,10,.35,.28],
	[.18,.01,.32,.17,.35,10,.12],
	[.13,.13,.29,.14,.28,.12,10]];

my $d4 = [
	[10,100,100,100,100,100,100],
	[.19,10,100,100,100,100,100],
	[.27,.31,10,100,100,100,100],
	[.08,.18,.26,10,100,100,100],
	[.33,.36,.41,.31,10,100,100],
	[.18,.01,.32,.17,.35,10,100],
	[.13,.13,.29,.14,.28,.12,10]];

#Single Linkage Matrix:http://people.revoledu.com/kardi/tutorial/Clustering/Numerical%20Example.htm
my $d5 = [
	[0.00,0.71,0.56,3.61,4.24,3.20],
	[0.71,0.00,4.95,2.92,3.54,2.50],
	[5.66,4.95,0.00,2.24,1.41,2.50],
	[3.61,2.92,2.24,0.00,1.00,0.50],
	[4.24,3.54,1.41,1.00,0.00,1.12],
	[3.20,2.50,2.50,0.50,1.12,0.00]];
#Complete Linkage Example:http://www.lx.it.pt/~afred/tutorials/B_Clustering_Algorithms.pdf pg7
my $d6 = [
	[1e10,4.00,11.7,20.0,21.5],
	[4.00,1e10,8.10,16.0,17.9],
	[11.7,8.10,1e10,9.80,9.80],
	[20.0,16.0,9.80,1e10,8.00],
	[21.5,17.9,9.80,8.00,1e10]];

# print "@{[fmd($D, [0,3,0])]}\n";
# ClusteringTree::print(upgma($D));
print "UPGMA:\n";
ClusteringTree::print(ClusteringAlg::upgma($d1));
print "\nWPGMA:\n";
ClusteringTree::print(ClusteringAlg::wpgma($d2));
print "\nSingle Linkage:\n";
ClusteringTree::print(ClusteringAlg::single_linkage($d5));
print "\nComplete Linkage:\n";
ClusteringTree::print(ClusteringAlg::complete_linkage($d6));
