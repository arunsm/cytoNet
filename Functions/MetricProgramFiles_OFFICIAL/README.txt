1. How to install the metrics code
2. How to run the metrics code
3. How to add a new metric



1. How to install the metrics code

To install the metrics code, move the Metrics folder to a machine with Matlab.
If you are going to run the metrics on the Davinci cluster, also copy the 
metrics.pbs file.


2. How to run the metrics code

  A. If you are NOT running on Davinci, do the following:
    1. Start Matlab and edit the frontend1.m file in the Metrics folder.
      a. Set the cellmasksuffix variable to be the common suffix for all 
         cell mask file names.
      b. Set the nucleusmasksuffix variable to be the common suffix for 
         all nucleus mask file names.
      c. Set the stainsuffix variable to be the common suffix for all stain
         mask file names.
      d. Set the stainorder variable to the stain order in the stain files.
    2. Run the fontend1 function passing it the following 3 arguments:
      a. As a string, the name of the single folder containing the mask and
         image files.  Each mask file contains a sequence of images, one per cell.
      b. As a string the name of the folder containing the individual metrics.
         Currently, this is 'IndivdualMetrics'.  If you wish to run just a single 
         metric, the second parameter should be the name of the metric file.
      c. As a string, the name to be given to the output file.

  B. If you ARE running on Davinci, do the following:
    1. Edit the frontend.m file making changes a-d described in section A above.
       You may find it easier to make the changes before transferring it to the
       Davinci cluster.
    2. Edit the metrics.pbs file and make the following changes:
      a. Refer to https://docs.rice.edu/confluence/display/ITDIY/Getting+Started+on+DAVinCI#GettingStartedonDAVinCI-PBSBatchScriptOptions
      b. Set the walltime.
      c. Set the queue name.
      d. Find the call to the Matlab addpath function and set its argument to 
         the name of the metrics folder.
      c. Find the call to the frontend function and set its 3 arguments as 
         described in section A above.
      d. Type the following at the Linux command line:

           qsub metrics.pbs

         Write down the job number returned by the qsub command
      e. When the job is completed, two files are created: metrics.pbs.eNNNN
         and metrics.pbs.oNNNN where NNNN is the job number.  Inspect these files
         to make sure that the metrics ran successfully.
      f. Learn how to get the status of a job and other details at:
         https://docs.rice.edu/confluence/display/ITDIY/Getting+Started+on+DAVinCI#GettingStartedonDAVinCI-SubmittingandMonitoringJobs

3. How to add a new metric
  Metrics functions take a single argument and return two values:

    function [header table] = mymetric(data)

  A. The data input variable is a Matlab struct.
     (http://www.mathworks.com/help/matlab/structures.html)
     The data variable contains the following fields:

    1. cellmaskdata - a cell array containing individual cell masks.
       (http://www.mathworks.com/help/matlab/matlab_prog/what-is-a-cell-array.html)
    2. nucleusmaskdata - a cell array containing individual nucleus masks.
    3. staindata - a cell array containing stain images.
    4. stainkey - a struct whose fieldnames are stain names containing the index of
       each stain in staindata

    Example of the beginning of a metric function:

      function [header table] = mymetric(data)
      % Get stain index of nuclei
      key = data.stainkey;
      nucleusindex = key.nucleus;

      % Get nuclei stain
      stains = data.staindata;
      nucleistain = stains(nucleusIndex);

      % Get cell masks
      cellmasks = data.cellmaskdata;
      numcells = numel(cellmasks);

      % Create header and table
      header = {'My Metric'};         <----- Note the curly braces
      table = zeros(numcells, 1);

      % Loop through each cell mask
      for i = 1:numel(cellmasks)
        onecellmask = cellmasks{i};   <----- Note the curly braces


  B. The metric function should place the metric values in the output variable
     table.  There should be a row for each cell.  If the function computes more
     than one metric, each different metric value is placed in a different column.

  C. The metric function should place the names of the individual metrics 
     (as strings) in the output variable header.  The header output variable should
     be a cell array.  Each column in the output variable table should have a name
     in the header output variable.

  Put the file containing the metrics function in the IndividualMetrics folder
  inside of the Metrics folder.  When the metrics are run, the metrics engine will
  find the file and use it.
