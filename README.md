# BridgePuzzle
Python program for solving a Bridge Puzzle (Hashiwokakero)

Files included:
  - BridgeObjects.py  -> The python file with the classes to solve a bridge puzzle
  - test.csv          -> A test puzzle as well as an example of how to fill out the puzzle
  - Output.csv        -> The solution to test.csv as solved by the program and confirmed by the source of the puzzle (Linkdoku iPhone app)

The file named BridgeObjects.py contains the three classes Node, Edge, BridgePuzzle that make up this type of puzzle. The BridgePuzzle class performs the main functionality

To solve one of these puzzles using this file, do the following:
- Create a csv file in the same directory as BridgeObjects.py
  - In this csv file, fill all locations with 0
  - Where there are nodes in the puzzle, place the number of connections that may be made to that node
- Import the BridgeObjects.py file
- Create a BridgePuzzle object with the file name of your previous csv as an input ( puzzle = BridgePuzzle("test.csv") )
  - Alternatively, you can do the following:
    - Create a BridgePuzzle object ( puzzle = BridgePuzzle() )
    - Call the read_file function ( puzzle.read_file("test.csv") )
    - Call the solve function ( puzzle.solve() ) where you may provide a file name to output the solved puzzle to. If no file name is specified, the solution will be written to Output.csv in the same directory as the python file
- The output will be written to a csv file in the same directory as the python file with the name Output.csv
  - To best read this file, I recommend centering the alignment in all cells and shrinking the cell sizes to just what is necessary to display the solved puzzle
