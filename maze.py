import numpy as np
import matplotlib.pyplot as plt
from random import choice, sample
from PIL.Image import open
# for creating a maze with r row and c column 
# an array[2r-1][2c-1] is needed in order to create walls
#
#
# an average performance up to 50 x 50 = 2500 pix
#
#
# creating a maze[ROW x COLUMN]

ROW = int(input('# Row    :\t'))
COL = int(input('# Column :\t'))
raw_maze = np.ones((2 * ROW - 1, 2 * COL - 1), dtype = int)

raw_maze[::2,::2] = 0

class runner:
    def __init__(self, maze):
        self.maze                       = maze
        
        self.previousRow                = -1
        self.previousColumn             = -1
        
        self.currentRow                 = 0
        self.currentColumn              = 0

        self.currentPosition            = [self.currentRow,self.currentColumn]
        self.previousPositions          = [self.currentPosition]
        self.transpasses                = []
        
    def go(self):
        getattr(self,choice(['changeRow', 'changeColumn']))()
    
    def changeRow(self):
        move                            = choice([-2, 2])
        if not [self.currentRow + move, self.currentColumn] in self.previousPositions and self.currentRow + move in range(0, len(self.maze)):
            self.previousRow            = self.currentRow
            self.currentRow             = self.currentRow + move
            self.transpass              = [[self.currentRow - move // 2, self.currentColumn]]
            self.currentPosition        = [[self.currentRow, self.currentColumn]]
            self.transpasses            += self.transpass
            self.previousPositions      += self.currentPosition
        """
        else:
            
            # to many repeating
            self.go()
            
            pass
        """
    def changeColumn(self):
        move                            = choice([-2, 2])
        if not [self.currentRow, self.currentColumn + move] in self.previousPositions and self.currentColumn + move in range(0, len(self.maze[0])):
            self.previousColumn         = self.currentColumn
            self.currentColumn          = self.currentColumn + move
            self.transpass              = [[self.currentRow, self.currentColumn - move // 2]]
            self.currentPosition        = [[self.currentRow, self.currentColumn]]
            self.transpasses            += self.transpass
            self.previousPositions      += self.currentPosition
        """
        else:
        
            # to many repeating
            self.go()
            
            pass
        """
    
x = runner(raw_maze)
jump = []
while len(x.previousPositions) != ROW * COL:
    if len(jump) > 4:
        if jump[0] == jump[1] and jump[1] == jump[2] and jump[2] == jump[3]:
            x.currentRow, x.currentColumn = choice(x.previousPositions)
        jump = []
    jump.append(len(x.previousPositions))
    x.go()


"""
#too slow to draw corridors
for r,c in x.previousPositions:
    raw_maze[r][c]  = 0
"""

for r,c in x.transpasses:
    raw_maze[r][c]  = 0

plt.imshow(raw_maze, cmap = plt.cm.binary)

# just a fancy naming
maze_name = '{} X {} - Maze ({}).jpg'.format(ROW, COL, ''.join(sample('cayicerdenizhoioynar',5)))
plt.imsave(maze_name, raw_maze,cmap       = plt.cm.binary)

# resizing the maze
maze                                      = open(maze_name)
maze.resize((ROW * 100, COL * 100)).save(maze_name)
plt.show()
