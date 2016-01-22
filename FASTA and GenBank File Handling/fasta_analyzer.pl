#!/usr/bin/perl
use strict;
use warnings;

#Purpose: Opens and analyzes FASTA format files
#Author: DAVID KWAN

my $help = <<HELP;
Opens and analyzes all FASTA format files in the subdirectory "data".
Possible command line arguments for this program are -h for help, -a for about, or -p to run the program.
HELP

my $lastmod = localtime ((stat $0)[9]);
my $about = <<ABOUT;
Author: David Kwan
Date Last Modified: $lastmod
ABOUT

if (scalar @ARGV != 1 || $ARGV[0] !~ /^(-h|-p|-a)$/) { #detects invalid arguments
    print "You have entered an invalid option.\n", $help; #shows help heredoc if argument is invalid
} elsif ($ARGV[0] eq '-h') {
    print $help; #show the help heredoc if argument is -h
} elsif ($ARGV[0] eq '-a') {
    print $about; #show the about heredoc if argument is -a
} elsif ($ARGV[0] eq '-p') { #run the program if argument is -p
    opendir my $data, "data" or die "Could not open directory 'data'\n";
    my @filenames = readdir $data; #read all files in the directory data
    closedir $data;

    my @fastafilenames = grep /\.fa/, @filenames; #search for FASTA format files
    my @sortedfastafilenames = sort @fastafilenames; #sort the FASTA format files alphanumerically

    foreach (@sortedfastafilenames) {
        open(my $fastafilehandle, "<" , "data/$_") or die "Could not open FASTA file\n"; #open each FASTA file one at a time
        my @lines = <$fastafilehandle>; #grab all the lines in the file
        close $fastafilehandle;

        chomp @lines;
        my $description = shift @lines; #isolate the description line
        my @descriptionlist = split /\|/, $description; #isolate each segment of the description
        my @accession_version = split /\./, $descriptionlist[3]; #isolate the accession number and the version number
        my $seqdata = join "", @lines; #combine all the sequence data into one continuous string

        my $counttotal = length $seqdata; #count the sequence length
        my $counta = ($seqdata =~ s/A/A/gi); #count the number of each base
        my $countc = ($seqdata =~ s/C/C/gi);
        my $countg = ($seqdata =~ s/G/G/gi);
        my $countt = ($seqdata =~ s/T/T/gi);

        print "*" x80, "\n"; #PRINT AN ANALYTICAL REPORT OF THE FASTA FILE
        print "filename: ", "$_", "\n";
        print "gi number: ", $descriptionlist[1], "\n";
        print "accession number: ", $accession_version[0], "\n";
        print "version: ", $accession_version[1], "\n";
        print "description: ", $descriptionlist[4], "\n";
        print "sequence length: ", $counttotal, "\n\n";
        printf ("%10s\t%10s\t%10s\n", "base", "count", "percentage");
        printf ("%10s\t%10s\t%5.2f\n", "a", "$counta", $counta/$counttotal*100);
        printf ("%10s\t%10s\t%5.2f\n", "c", "$countc", $countc/$counttotal*100);
        printf ("%10s\t%10s\t%5.2f\n", "g", "$countg", $countg/$counttotal*100);
        printf ("%10s\t%10s\t%5.2f\n", "t", "$countt", $countt/$counttotal*100);
    }

    print "*" x80, "\n"; #PRINT A SUMMARY OF THE PROGRAM EXECUTION
    print "Summary\n";
    print "Report printed on: ", my $localtime = localtime(), "\n";
    print "Program used: ", $0, "\n";
    print "Number of files read: ", scalar @filenames, "\n";
    print "Number of files processed: ", scalar @fastafilenames, "\n";
    print "Files processed: ", "\n";
    foreach (@sortedfastafilenames){
            print $_, "\n";
    }
    print "*" x80, "\n";
}
