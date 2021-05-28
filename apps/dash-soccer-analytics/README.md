# Soccer Match Analytics

## Overview

* `Soccer Match Analytics` helps coaches and analysts to analyze the events and view a digital recreation of action in a single match
  
## Usage
The app is actually comprised of two parts: a visual app and a game animation preparation script. This pre-preparation of animated game activity is necessary in order to speed up the graphical rendering process and minimize the amount of data processing downloading required to view a match. It is recommended to process no more than 25 minutes of a match at at time. Beyond this threshold it may be too difficult to create and render graphs.

Pre-processing of animated match data can be accomplished by doing the following:
- Executing the motion-graph.py script and selecting the time period of a match that you would like to pre-process (again in max 25 minute windows). You will need to select a .csv tracking file when executing the script.
- Save the resulting file in the data directory and name it using a .json file extension
- The file will now be visible and selectable within the app 
- The submit button must be selected to view the match

The event viewer is fairly self-explanatory and the user can select the team that they wish to see using the menus at the very top of the app. 

## Acknowledgements
With great gratitude I would like to thank Bruno Dadnino and his team at Metrica Sports (https://metrica-sports.com/about-us/) for the sample tracking and event data used in this application.
The original files are here: https://github.com/metrica-sports/sample-data. These files are part of Metrica's Elite Data product. So if you are subscribers to the Elite Data package you will be able to use those files with this application with very minimal changes required (they have been very lightly modified for the purposes of this app). Just copy the format used in the demo files (the changes/differences are very small). Better yet, do it programmatically and it's even easier!

