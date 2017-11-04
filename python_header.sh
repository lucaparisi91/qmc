for file in $(ls *.py)
do
    sed -i "1s:.*:#!/home/lucap/anaconda2/bin/python:" $file
done
