function name_chg_smi() {
  local fname=$1
  local line=$(sed -n '2 p' $fname)
  #local name=$(basename $fname)
  local smi=$(awk '{print $2}' <<<$line | cut -d":" -f2)
  local chg=$(awk '{print $1}' <<<$line | cut -d":" -f2)

  #printf "\"%s %s %s\"\n" $name $chg $smi
  printf "\"%s %s %s\"\n" $fname $chg $smi
}

function data_select() {
  local line=$(sed 's/\"//g' <<<$1)
  local lookfrom=$2
  local name=$(awk '{print $1}' <<<$line)
  local output=$(grep "$name" $lookfrom)
  if ! [ -z "$(echo $output)" ]; then
    printf "%s\n" "$1"
  fi

}

export -f name_chg_smi
export -f data_select

#XXXXXX="250227"

tot=name_chg_smi_neutral_tot.withRadical$XXXXXX

if [ -f $tot ]; then
  rm $tot
fi

parallel name_chg_smi :::: ./neutral_tot.withRadical$XXXXXX >$tot

lookfrom=./neutral_sel.woRadical$XXXXXX
sel=name_chg_smi_neutral_sel.woRadical$XXXXXX

parallel data_select {} $lookfrom :::: $tot >$sel
