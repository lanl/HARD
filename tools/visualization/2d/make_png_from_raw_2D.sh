#!/bin/bash

if [ $# -lt 4 ]
then
  echo
  echo 'Arguments missing.'
  echo '------------------'
  echo "Usage: $0 start_id end_id num_cells plot_column"
  echo "       start_id: (int) ID of first raw file to convert"
  echo "       end_id: (int) ID of last raw file to convert"
  echo "       num_cells: (int) number of cells in each direction (at the momemt: assumes same number in x and y direction)"
  echo "       plot column: (int) column in raw files to use for plotting value in z direction"
  echo "                    4: density"
  echo "                    5: pressure"
  echo "                    6: x component of velocity"
  echo "                    7: y component of velocity"
  echo "                    8: x component of momentum"
  echo "                    9: y component of momentum"
  echo "                   10: total energy"
  echo '  		   67: absolute value of velocity'
  echo '  		   89: absolute value of momentum'
  echo
  echo " Example: Plot density of 200 first time steps of a 128x128 hard run."
  echo "          $0 1 200 128 4" 
  echo
  exit 1;
fi

numcells=$3
column=$4
if [ $4 -eq 67 ]
then
  column='(sqrt(\$6*\$6+\$7*\$7))'
fi
if [ $4 -eq 89 ]
then
  column='(sqrt(\$8*\$8+\$9*\$9))'
fi

for i in `seq -f "%05g" $1 1 $2`
do
  if [ ! -f output-$i-2D-0.raw ]
  then
    echo "Can't find file output-$i-2D-0.raw. Abort."
    exit 1;
  fi
  echo -n "$i: time = "
  time=$(tail -n2 output-$i-2D-0.raw | head -n1 | awk '{print $1}')
  echo $time
  sed "s/XXXXX/$i/g; s/YYYYY/$numcells/g; s/ZZZZZ/$time/g; s/CCCCC/$column/g" output-XXXXX-2D-0.raw.gnuplot > output-current-2D-0.gnuplot
  gnuplot output-current-2D-0.gnuplot
done
