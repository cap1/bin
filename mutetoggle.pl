#!/usr/bin/env perl

#simple toggle for Speaker with amixer an console output

$totest = `/usr/bin/amixer get Speaker | /bin/grep ..Front.Left..Playback..off.`;
chomp($totest);
$compare= "  Front Left: Playback [off]";
if ( $totest eq $compare )
{
	system("/usr/bin/amixer -q set Speaker unmute");
	print(">> Speaker unmuted\n");
}
else
{
	system("/usr/bin/amixer -q set Speaker mute");
	print(">> Speaker muted\n");
}
