#!/bin/tcsh
set searchFor="Dashboards"
foreach myfile ($*)
  echo "File: $myfile"
  set myline=`otool -l $myfile |grep $searchFor`
  foreach line ($myline)
    if(`echo $line |grep -c $searchFor` > 0 ) then
      set libname=$line:t
      install_name_tool -change "$line" "@rpath/lib/Slicer-4.4/$libname" $myfile
    endif
  end
end
