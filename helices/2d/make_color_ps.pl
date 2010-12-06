#!/usr/bin/perl -w
use strict;
my @list = ('18S_rRNA.ps', '25S_rRNA_3p.ps', '25S_rRNA_5p.ps');
foreach my $postscript (@list) {
    my $input_text = $postscript;
    $input_text =~ s/ps$/txt/;
    Color_File($postscript);
}

sub Color_File {
    my $input_ps = shift;
    my $output_ps = $input_ps;
    $output_ps =~ s/\.ps/_colored\.ps/g;
    my $input_txt = $input_ps;
    $input_txt =~ s/ps$/txt/g;
    open(IN, "<${input_txt}") or print "Problem opening $input_txt $!\n";
    my @color_list = ();
    while (my $line = <IN>) {
	next if ($line =~ /^\#/);
	next if ($line =~ /^\s.$/);
	my ($num, $color) = split(/\s+/, $line);
	push(@color_list, $color);
    }
    close(IN);
    
    open(IN_PS, "<$input_ps");
    open(OUT_PS, ">$output_ps");
    my $count = undef;
  LINE: while (my $line = <IN_PS>) {
      chomp $line;
      print OUT_PS "$line\n";
      if ($input_ps eq "18S_rRNA.ps") {
	  if ($line eq "290.00 -105.33 290.00 -98.67 lwline") {
	      $count = 0;
	  } else {
	      next LINE if (!defined($count));	
	  }
      } elsif ($input_ps eq "25S_rRNA_3p.ps") {	
	  if ($line eq "-148.33 -1010.00 -141.67 -1010.00 lwline") {
	      $count = 0;
	  } else {
	      next LINE if (!defined($count));
	  }
      } elsif ($input_ps eq "25S_rRNA_5p.ps") {
	  if ($line eq "360.00 0.00 1.00 1.00 1.00 431.01 154.00 lwfarc") {
	      $count = 0;
	  } else {
	      next LINE if (!defined($count));
	  }
      }
      if (defined($count)) {
	  $count++;
      } else { ## $count is not defined
	  next LINE;
      }
      
      
      if ($count == 1) {
	  next LINE;
      } else {  ## count not 1 nor 0
	  my $list_counter = $count - 1;
	  my $color = $color_list[$list_counter];
	  next LINE if (!defined($color));
	  if (defined($color)) {
	      if ($color eq '0') {
		  ## 0 is black
		  print OUT_PS "0 0 0 setrgbcolor\n";
		  next LINE;
	      } elsif ($color eq '10') {
		  ## 10 is gray
		  print OUT_PS "0.6 0.6 0.6 setrgbcolor\n";
		  next LINE;
	      } elsif ($color eq '11') {
		  ## 11 is neon pink
		  print OUT_PS "0.85 0.30 0.64 setrbgcolor\n";
	      } elsif ($color eq '-4') {
		  ## -4 is purple
		  print OUT_PS "0.36 0.18 0.64 setrbgcolor\n";
	      } elsif ($color eq '-3') {
		  ## -3 is blue
		  print OUT_PS "0.08 0.25 1.0 setrbgcolor\n";
	      } elsif ($color eq '-2') {
		  ## -2 is greenblue
		  print OUT_PS "0.25 0.90 0.92 setrbgcolor\n";
	      } elsif ($color eq '-1') {
		  ## -1 is green
		  print OUT_PS "0.1 0.90 0.1 setrbgcolor\n";
	      } elsif ($color eq '1') {
		  ## 1 is yellow
		  print OUT_PS "0.9 0.9 0.15 setrbgcolor\n";
	      } elsif ($color eq '2') {
		  ## 2 is yelloworange
		  print OUT_PS "0.90 0.60 0.10 setrbgcolor\n";
	      } elsif ($color eq '3') {
		  ## 3 is orangered
		  print OUT_PS "0.92 0.34 0.08 setrbgcolor\n";
	      } elsif ($color eq '4') {
		  ## 4 is red
		  print OUT_PS "0.92 0.10 0.10 setrbgcolor\n";
	      } else {
		  ## default is... randomly chosen by me
		  print OUT_PS "0.57 0.08 0.32 setrbgcolor\n";
	      }
	  }
      }
  }
    close(IN_PS);
    close(OUT_PS);
}
