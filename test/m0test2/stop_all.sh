mkdir  backup
find . -mindepth 1 -maxdepth 1 -type d  -exec ./stop.sh {}  \;
