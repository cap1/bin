#!/bin/sh

status=`/usr/bin/xrandr --verbose | /bin/grep LVDS | /bin/sed "s/^[^ ]* [^ ]* [^ ]* ([^(]*) \([a-z]*\).*/\1/"`

if [ $status == "normal" ];
then
	/usr/bin/xrandr -o right
    /usr/bin/xsetwacom set "stylus" Rotate 1
    /usr/bin/xsetwacom set "cursor" Rotate 1
	/usr/bin/xsetwacom set "eraser" Rotate 1
else
	/usr/bin/xrandr -o normal
	/usr/bin/xsetwacom set "stylus" Rotate 0
	/usr/bin/xsetwacom set "cursor" Rotate 0
	/usr/bin/xsetwacom set "eraser" Rotate 0
fi
