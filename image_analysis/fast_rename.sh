for FILE in `find $1 -name '*.tif'`
do
#echo $FILE | sed 's/\S//g'

NEWNAME=`echo "$FILE" | sed -E 's/(\S+)\/.*$/schnuff/g'`
echo $NEWNAME;
done