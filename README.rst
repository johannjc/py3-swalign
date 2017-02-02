SWAlign
-------

Yet another smith waterson alignment code.

I copied the algorithm from: http://www.codesofmylife.com/2011/05/13/smith-waterman-algorithm-for-local-alignment-in-python/

And then tried to make it faster. 

At least on my laptop I was able to do about 300, ~160 bp alignments per second.

Oh yes, it also only works correctly on python3.


-------
This fork returns the maxScore as well - and fixed an issue with the setup.py to include the numpy include dirs
