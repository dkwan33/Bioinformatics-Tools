#!/usr/bin/perl
use strict;
use warnings;

#Purpose: Convert NCBI genbank (.gbk) files to FASTA (.fa) format for a given gene of interest
#Author: DAVID KWAN

my $help = <<HELP;
Convert all NCBI genbank (.gbk) files to FASTA (.fa) format in subdirectory "data" for a given gene of interest.
The gene must be given as a command line argument.
HELP

my $lastmod = localtime ((stat $0)[9]);
my $about = <<ABOUT;
Author: David Kwan
Date Last Modified: $lastmod
ABOUT

if (@ARGV == 1) { #if 1 argument is given, open the data directory, otherwise end program
    opendir my $data, "data" or die "Could not open directory 'data'\n";
    my @files = readdir $data; #read all files in the directory data

    if (my @gbkfiles = grep /\.gbk/, @files) { #if at least one .gbk file exists continue, otherwise end program
        my $filesuccess = 0;
        my @organism;
        my @seqaa;
        my @seqdna;

        foreach (@gbkfiles) { #for every .gbk file
            open(my $file, "<", "data/$_") or die "Could not open .gbk file\n"; #open the file
            my @lines = <$file>; #grab the lines from the file
            close $file; #close the file

            if ((grep /gene="$ARGV[0]"/, @lines) == 0) { #if gene isnt found, proceed to next file
                print "Could not find gene $ARGV[0] in file $_\n";
                next;
            }
            my $linesjoined = join "", @lines;
            $linesjoined =~ /.*organism="(.*?)".*?gene="$ARGV[0]".*?CDS.+?(complement)?\(?(\d+)\.\.(\d+).*?translation="(.*?)".*ORIGIN(.*)\/\//s; #regex that finds all useful information at once
            my $complement = 0;
            my ($organism, $index1, $index2, $seqaa, $seqdnaraw) = ($1, $3, $4, $5, $6); #use backreferences to pull out all useful data into usable variables
            if ($2) { #if the complement indicator is found
                $complement = 1; #set $complement to true
            }

            $organism =~ s/\s/_/g; #remove spaces in organism name
            $seqaa =~ s/\s//g; #remove spaces in AA sequence
            $seqdnaraw =~ s/\s|\d//g; #remove spaces and numbers in dna sequence

            my $seqdna = substr $seqdnaraw, $index1-1, $index2-$index1+1; #find gene sequence in the raw dna data

            if ($complement) { #if the complement indicator variable ($complement) is set to true
                $seqdna = reverse $seqdna; #change the dna seq to be the reverse complement
                $seqdna =~ tr/acgt/tgca/;
            }

            (my $seqdna80 = $seqdna) =~ s/(.{1,80})/$1\n/g; #space dna seq into FASTA format
            (my $seqaa80 = $seqaa) =~ s/(.{1,80})/$1\n/g; #space amino acid seq into FASTA format

            push @organism, $organism; #store all ready-to-print information for later use when printing
            push @seqdna, $seqdna80;
            push @seqaa, $seqaa80;

            $filesuccess++; #increase counter of files successfully read and analyzed
        }
        if (scalar @gbkfiles == $filesuccess) { #write report of gene to FASTA files if all files contained gene of interest
            open(my $dnaoutput, ">", "data/dna_out_$ARGV[0].fa") or die "could not open dna file to write to.\n";
            open(my $aaoutput, ">", "data/aa_out_$ARGV[0].fa") or die "could not open AA file to write to.\n";
            print "Your gene DNA and amino acid sequences have been written to dna_out_$ARGV[0].fa and aa_out_$ARGV[0].fa respectively in the subdirectory 'data.'\n";

            foreach (0..($filesuccess-1)) { #write the gene sequence to each FASTA file for each .gbk file
                print $dnaoutput ">$organism[$_]\n", "$seqdna[$_]";
                print $aaoutput ">$organism[$_]\n", "$seqaa[$_]";
            }
            close $dnaoutput;
            close $aaoutput;
        } else {
            print "No data was saved because gene '$ARGV[0]' was not found in one or more of the .gbk files.\n";
        }
    } else {
        print "There were no .gbk files found in the directory 'data'.\n";
    }
    closedir $data;
} else {
   print "Please enter ONE gene of interest as a command line argument.\n";
}
