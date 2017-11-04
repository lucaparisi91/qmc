cd $1
kill $(cat qmc.log | cut  -d'=' -f2)
