mkdir  backup
find . -mindepth 1 -maxdepth 1 -type d -not -name backup -exec mv  {} backup \;
