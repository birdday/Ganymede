# Ganymede

Solution to Ganymede take home exercise by Brian A. Day.

## Installation Instructions

Clone the repository. In a new environment (recommended), run `pip install -e .` from this directory. Note, I had to update pip for a successful install of the packages. Confirm successful install by running tests with `pytest`. Some tests throw a warning becuase of the packages used, but you should recieve no errors. Check out the examples.ipynb file for examples/solutions to the problems given below.

## Problem Statement

Of the many instruments and data that SSE’s come across, among the most prevalent are various forms of liquid chromatography (LC) machines. These instruments enable users to separate and precisely quantify constituents of mixtures across many different properties. In this take-home, you’ll write some code to parse and analyze some data output from an LC instrument.

1. Write a parser in Python that reads in the example data located [here](https://github.com/Tarskin/HappyTools/blob/master/Example%20Data/IgG%20Vtag%201_ACQUITY%20FLR%20ChA.txt) and creates a Python class that holds the following:

   - Metadata for the chromatogram run - injection information, chromatogram data information, and the signal parameter information section
   - Raw chromatogram data

2. Add a method to the class which determines peaks from the raw chromatogram data. This function should return the left threshold, right threshold, and height of the peak.
3. Add a method to the class which determines the elution volume represented by integrating each of the peaks.
4. Document your code, noting any caveats and limitations to what you have written.
   As you’re working through this assignment, consider how your code would look if presented to a client. Beyond functionality, documentation and code quality will be important attributes we look for in a strong response.
