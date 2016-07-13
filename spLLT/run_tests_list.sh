for mat in "${matrices[@]}"
do

    matrix=$(echo $mat | awk -F'matrix:' '{print $2}' | cut -f1 -d' ')
    nb=$(echo $mat | awk -F'nb:' '{print $2}' | cut -f1 -d' ')
    ne=$(echo $mat | awk -F'ne:' '{print $2}' | cut -f1 -d' ')


done
