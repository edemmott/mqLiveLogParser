# mqLiveLogParser
Matlab script for interpreting Maxquant.Live log files generated through use of the targeting app.

**Todo: parser is functional but slow (3-4 minutes): the final loop for extracting fields to table format should be converted to regexp.**

## How to call:
```matlab
outputTable = mqlLogParser('filename.txt');
```

The output table rows correspond to each time a row was observed in a survey scan. These can be used to investigate and further optimise targeting parameters, topN through cycle, as well as mzDeciviation, Intensity, RT and their correction factors throughout the run.

![example image](/img/example_output_plots.png)

Ed Emmott, May 2019, Northeastern U.

email: e.emmott@northeastern.edu
www: http://edemmott.co.uk
twitter: http://twitter.com/edemmott
