#old_file_name=K40.mac
new_file_name=test.mac

#old_file=OneKPlate
#new_file=LabFloor_extended
#old_file=0.463
#new_file=0.557
old_file=100000
new_file=10000000


#cat $old_file_name | sed 's/'"$old_file"'/'"$new_file"'/g' > $new_file_name;
#rm -r $old_file_name;
#mv $new_file_name $old_file_name

declare -a StringArray=("Co60.mac" "K40.mac" "NeutronTh232.mac" "NeutronU238.mac" "Pb_210_206.mac" "Th232.mac" "U238.mac" "U_238_210.mac")
#declare -a StringArray=("Co60.mac" "K40.mac" "NeutronTh232.mac")

# Iterate the string array using for loop
for old_file_name in ${StringArray[@]}; do
  cat $old_file_name | sed 's/'"$old_file"'/'"$new_file"'/g' > $new_file_name;
  rm -r $old_file_name;
  mv $new_file_name $old_file_name
done
