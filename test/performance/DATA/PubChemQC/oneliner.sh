XXXXXX="250227"
find "$(pwd -P)/neutral" -type f > "neutral_tot.withRadical$XXXXXX"
python ./radical_check.py > "neutral_tot.woRadical$XXXXXX"
python ./data_selector.py
source name_chg_smi.sh
