WD="/media/matvey/57a42bb6-576e-435f-a711-e6cb2026e82d/LshCas13a_RNA_cleavage/"
cd $WD
./LshCas13a_C3000/Scripts/raw_data_processing.sh
./LshCas13a_d10LVM/Scripts/raw_data_processing.sh
./LshCas13a_in_vitro_tRNAs/Scripts/raw_data_processing.sh
./LshCas13a_in_vitro_total_RNA/Scripts/raw_data_processing.sh

./LshCas13a_C3000/Scripts/read_mapping.sh
./LshCas13a_d10LVM/Scripts/read_mapping.sh
./LshCas13a_in_vitro_tRNAs/Scripts/read_mapping.sh
./LshCas13a_in_vitro_total_RNA/Scripts/read_mapping.sh

./LshCas13a_C3000/Scripts/return_fragment_coords_table.py &
./LshCas13a_d10LVM/Scripts/return_fragment_coords_table.py &
./LshCas13a_in_vitro_tRNAs/Scripts/return_fragment_coords_table.py &
./LshCas13a_in_vitro_total_RNA/Scripts/return_fragment_coords_table.py &
wait

./LshCas13a_C3000/Scripts/merge_ends_count_tables.py 
./LshCas13a_d10LVM/Scripts/merge_ends_count_tables.py
./LshCas13a_in_vitro_tRNAs/Scripts/merge_ends_count_tables.py 
./LshCas13a_in_vitro_total_RNA/Scripts/merge_ends_count_tables.py 

