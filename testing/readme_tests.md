# test folder for matmodfit 
All tests files should be located in this folder or its subdirectories. 
No other files or folders than those described in this readme should exist in this directory!

run_test_cases.py: Runs the test cases and makes a test_report.txt in the top directory. 
clean_test.py: Cleans the generated test results.
test_cases.txt: Contains description of each test case, used by run_test_cases.py

The report file can be renamed (see below) and moved to the TestReports folder if it should be saved.
It should never be committed in the version control without being moved and renamed!

## TestCaseFolder
-   Contains all the test input files and experiment files. Experiment files should start with the same name as the input file
-   Contains additional files (e.g. abaqus input files), these should also start with the same name as the matmodfit input file
-   Each test input should be named test_NN.inp where NN is the test number.

## TestReports
-   Contains generated test reports. 
-   Report files should be named: test_report_mmf_vX_Y_yyyy_mm_dd 
    * X_Y is the major and minor version number
    * yyyy_mm_dd is the date the test was run
    
## AbaqusTestGenerationScripts
-   Contains scripts to generate input files and experiment files which can be used to verify that matmodfit gives the same results as Abaqus