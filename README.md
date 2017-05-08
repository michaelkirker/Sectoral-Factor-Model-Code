# Sectoral-Factor-Model
Model code for the model described the paper "What Drives Core Inflation? A Dynamic Factor Analysis of Tradable and Nontradable Prices". RBNZ Working Paper DP2010/13

Rather than include all the code to handle both New Zealand and Australian data, and all the different graphs and tables found in the Working Paper, I have opted to simplify down the code in this repository to the code required to estimate the model, compute core inflation, and produce a few of the key graphs. This should hopefully be more accessible to the user.

## Disclaimer(s) ##

Please note that the model code presented here reflects an updated version of that discussed in the RBNZ Working paper DP2010/13. A number of small corrections were made to the code.

While the structure of the model in the code is the same as that used by the RBNZ for their Sectoral Factor Model, the code here does not represent the exact same code being used by the RBNZ, and so results produced by this code are likely to differ from those published by the RBNZ. 


## Contact details ##

Created by: Michael Kirker

Email: <mkirker@uchicago.edu>

Website: [http://michaelkirker.net](http://michaelkirker.net "http://michaelkirker.net")

Git repository: [https://github.com/michaelkirker/Sectoral-Factor-Model-Code]("https://github.com/michaelkirker/Sectoral-Factor-Model-Code")


## Repository structure ##

* /code/
	* Folder containing functions required to run the model.
* /documentation/ 
	* Collection of documentation providing more background on the code
* /input/
	* Folder containing the data input required to estimate and run the model
* /output/ 
	* Folder containing all the finalized output of the model.
* /temp/ 
	* Folder containing temporary files created by the model.
* batch.m
	* Batch file that executes the model




## How to use this code ##

The entire model can be run using the `batch.m` file. Each part of this batch file can be run independently if you have already run the previous sections. This is useful if you want to make changes to the code the analyses the output without having to re-estimate the model. Simply comment out the sections you do not wish to run.


The user can change settings related to the data etc inside of `batch.m`. Changes to the model structure (i.e. number of lags) can be made within the file `DFMcore.m`. User input should not be required anywhere else.

## Version History ##

* Current version (0.1)
	* This version of the model corrects a couple of coding typos in the working paper version of the model. These changes are yet to be fully documented here yet, however I will try to document these changes at some point in the future. The only noticeable difference is that the core inflation measure produced by this code is smoother than that in the Working Paper version (which had a typo when updating a variance term in the model).
	* To simplify things for the user, not all output from the paper is produced.
	* The model file (`DFMcore.m`) still needs cleaning up and further commenting. See the documentation folder for more information on the mathematical steps being taken in estimating the model.




## To-do ##

The following is a list of forthcoming changes

* Improve documentation within the code
	* Specific focus on detailing all the elements within each output structure
	* Add in code for all the tables and graphs in the paper.
	