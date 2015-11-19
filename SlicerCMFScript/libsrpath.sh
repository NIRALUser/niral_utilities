searchFor="Dashboards"
for myfile in $(ls *.dylib);
do
  echo "File: $myfile"
  for line in $(otool -l $myfile | grep $searchFor );
  do
    if [ `echo $line | grep -c $searchFor` -eq 1 ];
    then
      libname=$(basename $line)
      echo install_name_tool -change $line @rpath/lib/Slicer-4.4/$libname $myfile
      install_name_tool -change $line @rpath/lib/Slicer-4.4/$libname $myfile
    fi
  done
done