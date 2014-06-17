#!/usr/bin/perl -w
use strict;

# $Id: fasta2entropy.pl,v 1.2 2001/12/11 00:36:57 bruns Exp $
# $Header: /usr/data/cvs/sequoia1/fasta2entropy.pl,v 1.2 2001/12/11 00:36:57 bruns Exp $
# $Log: fasta2entropy.pl,v $
# Revision 1.2  2001/12/11 00:36:57  bruns
# Added command line argument to specify which sequence is the reference
# Other polishing of code and comments
#
# Revision 1.1  2001/12/11 00:06:45  bruns
# Initial crude script to compute position specific entropy/ID in alignments
#

# Read in a fasta file containing a multiple sequence alignment
# Output one position per line, with entropy or percent ID at each position

my $usage = "fasta2entropy.pl [reference sequence number] < family.fsa > entropy.table";

# Get command line arguments
# For now, just the reference sequence for ID computation
my $ref_seq_num = shift;
$ref_seq_num = 1 unless defined $ref_seq_num;
die unless $ref_seq_num >= 1;
$ref_seq_num --; # convert from 1-based array to 0-based array

# Storage arrays for the sequence entry parts
my @descriptions = ();
my @sequences = ();
my @ids = ();

# 1) Read in one sequence at a time
{
  local $/ = "\n>"; # split input on sequences
  while (<>) {
    my $fasta = $_; # Got one sequence entry

    # Get rid of ">", to put first and last entries on same footing with the rest
    chomp $fasta; # remove trailing "\n>"
    $fasta =~ s/^\s*>\s*//; # remove any leading ">" and associated spaces
    $fasta =~ s/\s*>\s*$//; # remove any trailing ">" and associated spaces

    # Separate the sequence from the header
    unless ($fasta =~ m/^
                        \s*(\S+) # non-space ID
                        \s*(.*)\n # rest of header line
                        ((.*\n)*.*) # sequence
                        $/x) { 
      die $fasta;
    }
    my $id = $1;
    my $description = $2;
    my $sequence = $3;
    $sequence =~ s/\s//g; # remove spaces, returns
    $sequence =~ tr/a-z/A-Z/; # convert to upper-case

    # Save this sequence into the pile
    push @sequences, $sequence;
    push @ids, $id;
    push @descriptions, $description;
  }
}

# Make sure there are at least two sequences
die unless scalar @sequences >= 2;
die unless scalar @sequences > $ref_seq_num;
# Verify that the sequence lengths are all the same
my $sequence_length = length($sequences[0]);
foreach my $sequence (@sequences) {
  die unless length($sequence) == $sequence_length;
}

my @columns = ();

# 2) Put the characters in columns
# simultaneously record background frequencies of residues
my %background = ();
my $total_chars = 0;
for (my $residue = 0; $residue < $sequence_length; $residue ++) {
  my $column = "";
  foreach my $sequence (@sequences) {
    my $char = substr($sequence, $residue, 1);
    if ($char =~ /^[A-Z]$/) { # only letters count for statistics
      $background{$char} ++;
      $total_chars ++;
    }
    $column .= $char;
  }
  push @columns, $column;
}

# 3) convert background counts to frequencies
foreach my $char (keys %background) {
  $background{$char} /= $total_chars;
}

# 3) Print out the results
foreach my $column (@columns) {
  print $column;

  # Entropy
  print "\t"; # tab delimited
  printf "%8.3f", -get_entropy(COLUMN => $column);

  # Percent identity
  print "\t"; # tab delimited
  printf "%8.3f", get_identity(COLUMN => $column, REFERENCE => $ref_seq_num);

  print "\n";
}

################################################################
# End of running program, the rest of this file is subroutines #
################################################################

# Subroutine to calculate the entropy of a string
sub get_entropy {
  my (%args) = @_; # read arguments into %args hash

  # string containing characters observed in the column
  die unless exists $args{COLUMN};
  die unless defined $args{COLUMN};
  my $column = $args{COLUMN};

  # Count the frequency of each character in the column
  my %char_counts = ();
  foreach my $char (split //, $column) {
    $char_counts{$char} += 1.0;
  }

  # Compute the entropy
  my $entropy = 0.0;
  foreach my $char (keys %char_counts) {
    next unless exists $background{$char}; # Ignore characters ignored for background

    # convert raw count to probability
    my $counts = $char_counts{$char};
    my $P = ($counts - 1) * (1.0 / (length($column) - 1)); # Observed probability

    next if $P <= 0;
    my $Q = $background{$char}; # background probablity
    $entropy -= ($P * log($P/$Q));
  }

  # convert to units of "bits"
  $entropy /= log(2.0);

  return $entropy;
}

# Subroutine to get percent identity of a string, with respect to one position
# (properly, this should be without considering self identity)
sub get_identity {
  my (%args) = @_; # read arguments into %args hash

  # string containing characters observed in the column
  die unless exists $args{COLUMN};
  die unless defined $args{COLUMN};
  my $column = $args{COLUMN};

  # which sequence is the reference?
  die unless exists $args{REFERENCE};
  die unless defined $args{REFERENCE};
  my $reference = $args{REFERENCE};

  my $percent_id = 0.0;
  my $reference_char = substr($column, $reference, 1);
  foreach my $char (split //, $column) {
    next unless exists $background{$char}; # Ignore characters ignored for background
    if ($char eq $reference_char) {
      $percent_id ++;
    }
  }
  
  # convert raw counts to percent
  if ($percent_id > 0) {
    $percent_id = ($percent_id - 1) * 100.0 * 1.0 / (length($column) - 1);
  }

  return $percent_id;
}

