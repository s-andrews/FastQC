FastQC - A Quality Control application for FastQ files
------------------------------------------------------

Most high throughput sequencers generate output in FastQ format.  This
format combines the base calls for the sequence which was generated with
an encoded quality value for each base which says how confident the
sequencer was that the base call generated was correct.

Before proceeding with the analysis of a sequence data set it is
a good idea to do some basic quality control checks on the raw data
to ensure that there are no hidden problems which might be more
difficult to detect at a later stage.

FastQC is an application which takes a FastQ file and runs a series
of tests on it to generate a comprehensive QC report.  This will
tell you if there is anything unusual about your sequence.  Each
test is flagged as a pass, warning or fail depending on how far it
departs from what you'd expect from a normal large dataset with no
significant biases.  It's important to stress that warnings or even
failures do not necessarily mean that there is a problem with your
data, only that it is unusual.  It is possible that the biological
nature of your sample means that you would expect this particular
bias in your results.

FastQC can be run either as an interactive graphical application 
which allows you to view results for multiple files in a single
application.  Alternatively you can run the program in a non
interactive way (say as part of a pipeline) which will generate
an HTML report for each file you process.

FastQC is a cross-platform application, written in java.  In theory it
should run on any platform which has a suitable java runtime environment.
Having said that we've only tested in on Windows, MacOSX and Linux
running the Oracle v1.6 to 1.8 JREs.  Please let us know what happened if
you try running it on other platforms / JREs.  Please see the detailed
instructions in the INSTALL.txt document to tell you how to get a 
suitable java version to run FastQC on your system.

If you have any comments about FastQC we would like to hear them.  You
can either enter them into the github bug tracker at:

https://github.com/s-andrews/FastQC/issues/

..or send them directly to simon.andrews@babraham.ac.uk.
