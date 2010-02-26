#!/usr/bin/env perl
use strict;
use warnings;

#####
#small hack for xrandr for a second monitor
#####

my $param = "null";				
if ( $ARGV[0] ) {
	$param = shift(@ARGV);
}

my $status = `/usr/bin/xrandr`;
 
if ( $param eq "clone"and $status =~ /VGA1 connected/ ) {
     system("/usr/bin/xrandr --output LVDS1 --auto --output VGA1 --auto");
}
elsif ( $param eq "on"and $status =~ /VGA1 connected/ ) {
	system("/usr/bin/xrandr --output LVDS1 --auto --output VGA1 --auto --right-of LVDS1");
}
elsif ( $param eq "off" and $status =~ /VGA1 connected/ ) {
	system("/usr/bin/xrandr --output VGA1 --off");
}
else {
		print ">> somethings wrong here ...\n";
}
