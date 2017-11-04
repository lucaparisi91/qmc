find . -mindepth 1 -maxdepth 1 -type d -not -name backup -exec sh -c 'rm $0/*.dat' '{}' \;
