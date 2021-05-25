# Interactive-Visulizations-with-Plotly-Dash

# VISUALIZATION OF END GAME CHESS PIECES

## Motivation

Chess is one of the oldest games in history (reportedly invented in India, around the 7th centuryAD) that is still played today. As society evolved, so did the game, and how it is played. The advent of the Internet thrusted the game into an entirely new medium, enabling players to connect regardless of their location. Recent events like the global pandemic (which forced society to adapt their lifestyle) and the overwhelming popularity of the Netflix show “The Queen’s Gambit” revitalized the sport, as the number of played games online has skyrocketed in the past year. 

It is commonly said one can go their entire life without ever playing the same game of chess twice. However, with only 64 squares on the board, players will likely face similar situations very often. Their pieces too, as their movement patterns are defined by specific rules. Considering such a wide range of possible piece journeys, is there common ground at the end of that road? For this project, we have decided to answer such questions and visualize the most common position for the pieces in the board at the end of a game, via a heatmap. Does it change much depending on the player’s skill level? How much do pieces travel, considering games with different move length? How do these hotspots change for either player for different outcomes? For example, is there a square in which the King is more often checkmated?

## Dataset Description

The dataset was produced from data obtained from the Lichess open database [1], which hosts millions of games in PGN format, which have been parsed and stored as a .csv file. Due to size constraints and limitations, the PGN corresponding the month of April 2017, was chosen, and a sample of 5000 games were extracted from it. For this purpose, the pandas library was used in conjunction with the python-chess library. Information derived from the game includes winners, payer’s elo rating, played moves, chess piece positions, time control, and game types. 

## Visualization and Interactive Choices

Design wise, the presentation is meant to be simple but intuitive: a single chessboard where hotspots can be visualized was construed as the main piece of the visualization, while a range of interactable components are employed to add interactivity and functionality. Link to the app: https://chess-checkmate-king.herokuapp.com

## Technical Aspects

Accessing the dataset was the first challenge that had to be overcome, because of the large file size. Therefore, analyzing data over time was out of the question. Post extraction, the next challenge is in accessing, parsing, and extracting relevant information from the respective files. Chess games are typically encoded in Portable Games Notation (PGN) files which require programmatic APIs and libraries such as Python-chess to draw information from it and transform it into a plottable-ready structure. For this purpose, code was adapted from similar projects, and created from scratch, so that a new csv file that could be read and easily transformed was produced. 

Deploying an app online through Heroku is a herculean task that forces additional constraints on the layout and available options regarding interaction and visualization. Specifically, the memory requirements for web deployment required us to significantly downscale the dataset, choosing between additional filters/interactivity and sample size, because it was generally running out of memory during deployment. Making the visualizations smooth and fluid, as well as animated, was as challenging as it was rewarding. Specific techniques had to be employed under the hood, since plotly does not natively support animations for heatmaps which we had initially employed. Therefore, while the result technically looks like a chessboard, underneath it is a combination of three different scatterplots, with carefully manipulated axes and value outputs.

During our testing of the deployed web app, we were faced with issues related to updates to certain filters (which appeared to change randomly to different sets of choices), that we could not reproduce locally, or even a short while after redeployment. One possible hypothesis is that our input and output pipeline in Dash is not efficient, and the assignment of some variables fails in this process. Regardless, the output is still accurate. The app will also run into errors if there are no valid games to display in the stacked bar after filtering, which is a possibility if a user decides to seek out outliers, such as games with a very high number of moves. Part of this is related to how plotly updates and displays such graphs. 

For this reason, and because the interest of the visualization is on viewing a vast quantity of games, a catch-all feature for this issue was not prioritized.

## Discussion

We have deployed an app that allows us to visualize heatmap of pieces at the end of the game, as well as filter for specific conditions, game states, and pieces. Besides the dataset size limitation mentioned previously, the visualization is limited by the number of pieces’ heatmaps one can visualize at once, as well as being confined to end of game positions. Another limitation pertains the automatic scaling of elements in specific screens, which may break the order and shape of the visualization. In the future, one can expand the depth of this visualization by adding an additional chessboard so that both white and black pieces can be visualized at once (or other combination of pieces). 
One can also analyze a path for a specific piece for the entire game, or dynamically change the visualization for a specific point in the game. These two expansions would require the extraction and curated parsing of substantially more information from PGN’s. However, some of these visualizations already exist in the form of some published projects, such as [2] and [3].

## References

Go here.
[1] The Lichess Open Database. 2021 [online] Available at
<https://database.lichess.org/>[Accessed 5 April 2021].
[2] Kaggle.com. 2021. Chess games analysis. [online] Available at:
<https://www.kaggle.com/jac08h/chess-games-analysis> [Accessed 5 April 2021].
[3] ebemunk. 2016. A Visual Look at 2 Million Chess Games. Available at:
<https://blog.ebemunk.com/a-visual-look-at-2-million-chess-games/>[Accessed 5 April 2021].

