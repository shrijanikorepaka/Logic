from pysat.formula import CNF
from pysat.solvers import Solver
from typing import List

def solve_sudoku(grid: List[List[int]]) -> List[List[int]]:
    n = len(grid)
    ln = int(n ** 0.5)
    cnf = CNF()

    def P(i, j, d):
        # variables are 1..n for digits
        return i * n * n + j * n + d

    # Row, Col, and Cell constraints in one loop
    for i in range(n):
        for j in range(n):
            row_clause = []
            col_clause = []
            cell_clause = []
            for k in range(1, n + 1):
                # at least one per row, column, and cell
                row_clause.append(P(i, k - 1, j + 1))   # row i has digit j+1
                col_clause.append(P(k - 1, i, j + 1))   # col i has digit j+1
                cell_clause.append(P(i, j, k))          # cell (i,j) some digit

                # pairwise "at most one"
                for k2 in range(k + 1, n + 1):
                    cnf.append([-P(i, j, k), -P(i, j, k2)])   # cell only 1 digit
                    cnf.append([-P(i, k - 1, j + 1), -P(i, k2 - 1, j + 1)]) # row unique
                    cnf.append([-P(k - 1, i, j + 1), -P(k2 - 1, i, j + 1)]) # col unique

            cnf.append(row_clause)
            cnf.append(col_clause)
            cnf.append(cell_clause)

    # Subgrid constraints
    for d in range(1, n + 1):
        for br in range(ln):
            for bc in range(ln):
                cells = []
                for r in range(br * ln, (br + 1) * ln):
                    for c in range(bc * ln, (bc + 1) * ln):
                        cells.append((r, c))
                cnf.append([P(r, c, d) for (r, c) in cells])   # at least once
                for i1 in range(len(cells)):
                    for i2 in range(i1 + 1, len(cells)):
                        r1, c1 = cells[i1]
                        r2, c2 = cells[i2]
                        cnf.append([-P(r1, c1, d), -P(r2, c2, d)]) # at most once

    # Encode givens
    for i in range(n):
        for j in range(n):
            if grid[i][j] != 0:
                cnf.append([P(i, j, grid[i][j])])

    # Solve
    with Solver(name="glucose3") as solver:
        solver.append_formula(cnf)
        if not solver.solve():
            raise ValueError("No solution exists for this Sudoku.")
        model = set(solver.get_model())

    # Decode
    solution = [[0] * n for _ in range(n)]
    for r in range(n):
        for c in range(n):
            for d in range(1, n + 1):
                if P(r, c, d) in model:
                    solution[r][c] = d
                    break

    return solution
