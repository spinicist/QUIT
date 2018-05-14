#
# Regular cron jobs for the quit package
#
0 4	* * *	root	[ -x /usr/bin/quit_maintenance ] && /usr/bin/quit_maintenance
