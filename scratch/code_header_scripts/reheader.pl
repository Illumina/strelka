#!/usr/bin/env perl

use strict;
use warnings;

# new_header file
#

scalar(@ARGV) == 1 || die("usage: $0 new_header < file");

my $newhead = $ARGV[0];

open(my $HFH,"< $newhead") || die("Can't open $newhead");
my $newstr = "";
$newstr .= $_ while(<$HFH>);
close($HFH);

my $isfirst=1;
my $count=0;
while(<STDIN>) {
  if($isfirst) {
    if(/^\/\//) {
      $count++;
      next;
    }
    die("No header to replace!") unless($count);
    print $newstr;
    $isfirst=0;

    # if first line after comments is empty, then skip it:
    next if(/^\s*$/);
  }
  print $_;
}

die("No data!") if($isfirst);

