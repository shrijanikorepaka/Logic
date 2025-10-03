"""
Sokoban Solver using SAT (Boilerplate)
--------------------------------------
Instructions:
- Implement encoding of Sokoban into CNF.
- Use PySAT to solve the CNF and extract moves.
- Ensure constraints for player movement, box pushes, and goal conditions.

Grid Encoding:
- 'P' = Player
- 'B' = Box
- 'G' = Goal
- '#' = Wall
- '.' = Empty space
"""

from pysat.formula import CNF
from pysat.solvers import Solver

# Directions for movement
DIRS = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}


class SokobanEncoder:
    def __init__(self, grid, T):
        """
        Initialize encoder with grid and time limit.

        Args:
            grid (list[list[str]]): Sokoban grid.
            T (int): Max number of steps allowed.
        """
        self.grid = grid
        self.T = T
        self.N = len(grid)
        self.M = len(grid[0])

        self.goals = []
        self.boxes = []
        self.player_start = None

        # TODO: Parse grid to fill self.goals, self.boxes, self.player_start
        self._parse_grid()
        self.var_map={}
        self.var_cnt=0
        self.num_boxes = len(self.boxes)
        self.cnf = CNF()

    def _parse_grid(self):
        """Parse grid to find player, boxes, and goals."""
        # TODO: Implement parsing logic
        for row in range(self.N):
            for col in range(self.M):
                if self.grid[row][col]=='P':
                    self.player_start=(row,col) #players starting pos
                elif self.grid[row][col]=='B': #add to the list of boxes
                    self.boxes+=[(row,col)]
                elif self.grid[row][col]=='G':# add to the list of goals
                    self.goals+=[(row,col)] 
                 

    # ---------------- Variable Encoding ----------------
    def var_player(self, x, y, t):
        """
        Variable ID for player at (x, y) at time t.
        """
        # TODO: Implement encoding scheme
        if not self.isFree(x,y):
            return None
        return self._var('P',x,y,t)
    
    
    def isFree(self,x,y):
        #to check whether the position is inside the grid.
        return (0<=x<self.N )and (0<=y<self.M) and (self.grid[x][y]!='#')
    
    def var_box(self, b, x, y, t):
        """
        Variable ID for box b at (x, y) at time t.
        """
        # TODO: Implement encoding scheme
        if not self.isFree(x,y):
            return None
        return self._var('B',b,x,y,t)

    def _var(self,kind,*args):
        #assign a unique integer(+ve).
        key =(kind,)+args
        if key not in self.var_map:
            self.var_cnt+=1
            self.var_map[key]=self.var_cnt
        return self.var_map[key]
    
    
    
    # ---------------- Encoding Logic ----------------
    def encode(self):
        """
        Build CNF constraints for Sokoban:
        - Initial state
        - Valid moves (player + box pushes)
        - Non-overlapping boxes
        - Goal condition at final timestep
        """
        # TODO: Add constraints for:
        # 1. Initial conditions
        #set players's intital pos
        for x in range(self.N):
            for y in range(self.M):
                p=self.var_player(x,y,0)
                if p:
                    if (x,y)==self.player_start:
                        self.cnf.append([p])
                    else:
                       self.cnf.append([-p])
        #set box's initial positon
        for b in range(self.num_boxes):
            for x in range(self.N):
                for y in range(self.M):
                    b0=self.var_box(b,x,y,0)
                    if b0:
                        if (x,y)==self.boxes[b]:
                            self.cnf.append([b0])
                        else:
                            self.cnf.append([-b0])
       # 2. Player movement
        for t in range(self.T+1):
            player_vars=[]
            #player is in exactly one place
            for x in range(self.N):
                for y in range(self.M):
                    if self.isFree(x,y) and self.var_player(x,y,t):
                        player_vars+=[self.var_player(x,y,t)]
            self.cnf.append(player_vars) #atleast one
            for i in range(len(player_vars)): 
                for j in range(i+1,len(player_vars)):
                    self.cnf.append([-player_vars[i],-player_vars[j]])#atmost one
        
            #Each box is in exactly one place
            for b in range(self.num_boxes):
                box_vars=[]
                for x in range(self.N):
                    for y in range(self.M) :
                        if self.isFree(x,y) and self.var_box(b,x,y,t):
                            box_vars+=[self.var_box(b,x,y,t)]
            self.cnf.append(box_vars) #Atleast one
            for i in range(len(box_vars)):
                for j in range(i+1,len(box_vars)):
                    self.cnf.append([-box_vars[i],-box_vars[j]]) #Atmost one

        # 3. Box movement (push rules)
          #A box stays unless pushed 
        for t in range(self.T):
            for b in range(self.num_boxes):
                for x in range(self.N):
                    for y in range(self.M):
                        if self.isFree(x, y):
                            b_t = self.var_box(b, x, y, t)
                            b_t1 = self.var_box(b, x, y, t + 1)
                            if b_t and b_t1:
                                push_condn = [] #possible  push condtn
                                for dx, dy in DIRS.values():
                                    px, py = x - dx, y - dy
                                    if self.isFree(px, py):
                                        p_behind = self.var_player(px, py, t)
                                        p_on = self.var_player(x, y, t + 1)
                                        if p_behind and p_on: #player is behind the box & moves into the it
                                            mv = self._var('push', b, x, y, t, px, py)
                                            self.cnf.append([-mv, p_behind]) 
                                            self.cnf.append([-mv, p_on])
                                            self.cnf.append([mv, -p_behind, -p_on])
                                            push_condn.append(mv)
                                self.cnf.append([-b_t, b_t1] + push_condn)
            #Player moves and box pushing
            for x in range(self.N):
                for y in range(self.M):
                    if self.isFree(x, y):
                        p_t = self.var_player(x, y, t)
                        if p_t:
                            possible_moves = [] #check possible player moves
                            
                            for dx, dy in DIRS.values():
                                nx, ny = x + dx, y + dy
                                if self.isFree(nx, ny):
                                    p_t1 = self.var_player(nx, ny, t + 1)
                                    if p_t1:
                                        possible_moves += [p_t1]
                            self.cnf.append([-p_t] + possible_moves)
                            
                            for dx, dy in DIRS.values():
                                nx, ny = x + dx, y + dy
                                if not self.isFree(nx, ny):
                                    continue
                                p_t1 = self.var_player(nx, ny, t + 1)
                                if not p_t1:
                                    continue
                                bx, by = nx + dx, ny + dy #pos where box would be pushed
                                blocked_by_wall = not self.isFree(bx, by)
                                for b1 in range(self.num_boxes):
                                    b_in_way = self.var_box(b1, nx, ny, t)
                                    if blocked_by_wall:
                                        if b_in_way:
                                            self.cnf.append([-p_t, -b_in_way, -p_t1]) #cant move if there is a wall
                                    else:
                                        for b2 in range(self.num_boxes):
                                            b_at_dest = self.var_box(b2, bx, by, t)
                                            if b_in_way and b_at_dest: #check if it is free 
                                                self.cnf.append([-p_t, -b_in_way, -b_at_dest, -p_t1])
                                 #change box pos
                                if not blocked_by_wall:
                                    for b in range(self.num_boxes):
                                        b_t = self.var_box(b, nx, ny, t)
                                        b_t1_new = self.var_box(b, bx, by, t + 1)
                                        if b_t and b_t1_new:
                                            self.cnf.append([-p_t, -p_t1, -b_t, b_t1_new])
                                                                                                       
        # 4. Non-overlap constraints
            #No 2  obj can be in the same cell 
        for t in range(self.T+1):              
            for x in range(self.N):
                for y in range(self.M):
                    if not self.isFree(x,y):
                        continue
                    #NO overlap b/w player and box
                    p_var=self.var_player(x,y,t)
                    if p_var:
                        for b_ind in range(self.num_boxes):
                            b_var=self.var_box(b_ind,x,y,t)
                            if b_var:
                                self.cnf.append([-p_var,-b_var])
                        #No overlap b/w boxes
                    for b1 in range(self.num_boxes):
                        for b2 in range(b1+1,self.num_boxes):
                            b1_var=self.var_box(b1,x,y,t)
                            b2_var=self.var_box(b2,x,y,t)
                            if b1_var and b2_var:
                                self.cnf.append([-b1_var,-b2_var])
                            
        # 5. Goal conditions
        for b in range(self.num_boxes):
            goal_vars = []
            for gx, gy in self.goals:
                if self.isFree(gx, gy):
                    g_var = self.var_box(b, gx, gy, self.T)
                    if g_var:
                        goal_vars += [g_var]
            self.cnf.append(goal_vars)
        # 6. Other conditions
        pass
        return self.cnf


def decode(model, encoder):
    """
    Decode SAT model into list of moves ('U', 'D', 'L', 'R').

    Args:
        model (list[int]): Satisfying assignment from SAT solver.
        encoder (SokobanEncoder): Encoder object with grid info.

    Returns:
        list[str]: Sequence of moves.
    """
    N, M, T = encoder.N, encoder.M, encoder.T

    # TODO: Map player positions at each timestep to movement directions
    if not model:
        return -1
    rev_map={a:b for b,a in encoder.var_map.items()}
    player_pos=[None]*(encoder.T+1)
    for var in model:
        if var>0 and var in rev_map:
            kind,*args=rev_map[var]
            if kind=='P':
                x,y,t=args
                player_pos[t]=(x,y)
    moves=[]
    for t in range(encoder.T):
        if player_pos[t] is None or player_pos[t+1] is None:
            return -1
        x1,y1=player_pos[t]
        x2,y2=player_pos[t+1]
        dx,dy=x2-x1,y2-y1
        move_found=False
        for move_char,(mx,my) in DIRS.items():
            if(dx,dy)==(mx,my):
                moves.append(move_char)
                move_found=True
                break
        if not move_found:
            return -1
    return moves   



def solve_sokoban(grid, T):
    """
    DO NOT MODIFY THIS FUNCTION.

    Solve Sokoban using SAT encoding.

    Args:
        grid (list[list[str]]): Sokoban grid.
        T (int): Max number of steps allowed.

    Returns:
        list[str] or "unsat": Move sequence or unsatisfiable.
    """
    encoder = SokobanEncoder(grid, T)
    cnf = encoder.encode()

    with Solver(name='g3') as solver:
        solver.append_formula(cnf)
        if not solver.solve():
            return -1

        model = solver.get_model()
        if not model:
            return -1

        return decode(model, encoder)


