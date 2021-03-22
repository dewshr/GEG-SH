# Different case examples

## Running the program using only input file (case 1)
The program can be run using only input file using `python safe_harbor.py -i example_data/test_data.txt` or you can also pass output file name as `python safe_harbor.py -i example_data/test_data.txt -fname case1_results.csv`. In both cases output result will be same, only difference is the output file name. The final filtered result will be in file **filtered_case1_results.csv**, and **case1_results.csv** will contain unfiltered raw data informations.

The example file also contains FDR and AF columns, but the parameters are not passed while running the program. So, the program will not filter based on those columns and will generate a warining message as in picture below:

<img width="1071" alt="warning_message_safe_harbor" src="https://user-images.githubusercontent.com/22225447/112021513-5e667a80-8aff-11eb-95cd-b3c604659a15.png">

## Running the program passing FDR and allele frequency (AF) threshold values (case 2)
The user can provide FDR and AF thershold using `python safe_harbor.py -i example_data/test_data.txt -t .1 -af .1 -fname case2_results.csv`. This command will filter out the variants by removing the ones with fdr value and allele frequency less than given threshold.
