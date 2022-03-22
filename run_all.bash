./0_1_get_event_data.py
# 0.5_screen_event_data.py write bad-quality-stations.txt and rerun 0_1_get_event_data.py
#./ 0_1_get_event_data.py
cd data
../2_process_data_to_sac.py
bash fk_bash.cmd
cd ..
./3_prepare_for_gcap.py
bash cap_auto.bash
