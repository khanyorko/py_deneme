class matrix:
    def __init__(self, data):
        self.data   = [data] if type(data[0]) != list else data
        #if not all(list(map(lambda x: len(x) == len(self.data[0]), self.data))): raise ValueError("Matrix shape is not OK")
        self.row    = len(self.data)
        self.col    = len(self.data[0])
        self.shape  = (self.row, self.col)
#==================================================================================
    def __null__(self, shape, key = 1, diag = False):
        row, col = shape
        if row != col and diag: raise ValueError("Only square matrices could be diagonal")
        data = [[(key,(0, key)[i == j])[diag] for i in range(col)] for j in range(row)]
        return data
#==================================================================================
    def __str__(self):
        return "{}".format(self.data)
#==================================================================================
    def __repr__(self):
        str_ = ""
        for row in range(self.row):
            col_ = []
            for col in range(self.col):
                col_.append(str(self.data[row][col]))
            str_ += "\t".join(col_) + "\n"
        return "Matrix({}x{})\n".format(self.row, self.col) + str_
#==================================================================================
    def __add__(self, other):
        if self.shape != other.shape: raise ValueError("Matrix shapes are not matched")
        data = self.__null__(self.shape, key = 0)
        for row in range(self.row):
            for col in range(self.col):
                data[row][col] = self.data[row][col] + other.data[row][col]
        return matrix(data)
#==================================================================================
    def __sub__(self, other):
        if self.shape != other.shape: raise ValueError("Matrix shapes are not matched")
        data = self.__null__(self.shape, key = 0)
        for row in range(self.row):
            for col in range(self.col):
                data[row][col] = self.data[row][col] - other.data[row][col]
        return matrix(data)
#==================================================================================
    def __mul__(self, other):
        if type(other) != matrix: raise ValueError("Non-matrix type variables must be placed at the left side of the matrix \n\t\t------> const * A(mxk)")
        elif self.col != other.row: raise ValueError("Matrix shapes did not satisfy the criteria of C(mxn) = A(mxk) * B(kxn)")
        data = self.__null__((self.row, other.col), key = 0)
        for row in range(self.row):
            for col in range(other.col):
                for mul in range(self.col):
                    data[row][col] += self.data[row][mul] * other.data[mul][col]
        return matrix(data)
#==================================================================================
    def __rmul__(self, const):
        data = self.__null__(self.shape, key = 0)
        for row in range(self.row):
            for col in range(self.col):
                data[row][col] = self.data[row][col] * const
        return matrix(data)
#==================================================================================
    def __pow__(self, const):
        data = matrix(self.data)
        for i in range(const-1):
            data *= data
        return data
#==================================================================================
    def __pos__(self):
        data = self.__null__(self.shape, key = 0)
        for row in range(self.row):
            for col in range(self.col):
                data[row][col] = +self.data[row][col]
        return matrix(data)
#==================================================================================
    def __neg__(self):
        data = self.__null__(self.shape, key = 0)
        for row in range(self.row):
            for col in range(self.col):
                data[row][col] = -self.data[row][col]
        return matrix(data)
#==================================================================================
    def __abs__(self):
        data = self.__null__(self.shape, key = 0)
        for row in range(self.row):
            for col in range(self.col):
                data[row][col] = abs(self.data[row][col])
        return matrix(data)
#==================================================================================
    def __int__(self):
        data = self.__null__(self.shape, key = 0)
        for row in range(self.row):
            for col in range(self.col):
                data[row][col] = int(self.data[row][col])
        return matrix(data)
#==================================================================================
    def __float__(self):
        data = self.__null__(self.shape, key = 0)
        for row in range(self.row):
            for col in range(self.col):
                data[row][col] = float(self.data[row][col])
        return matrix(data)
#==================================================================================
    def __and__(self, other):
        data = self.__null__(self.shape, key = 0)
        for row in range(self.row):
            for col in range(self.col):
                data[row][col] = self.data[row][col] & other.data[row][col]
        return matrix(data)
#==================================================================================
    def __or__(self, other):
        data = self.__null__(self.shape, key = 0)
        for row in range(self.row):
            for col in range(self.col):
                data[row][col] = self.data[row][col] | other.data[row][col]
        return matrix(data)
#==================================================================================
    def __lt__(self, const):
        data = self.__null__(self.shape, key = 0)
        for row in range(self.row):
            for col in range(self.col):
                data[row][col] = self.data[row][col] < const
        return matrix(data)
#==================================================================================
    def __le__(self, const):
        data = self.__null__(self.shape, key = 0)
        for row in range(self.row):
            for col in range(self.col):
                data[row][col] = self.data[row][col] <= const
        return matrix(data)
#==================================================================================
    def __eq__(self, const):
        data = self.__null__(self.shape, key = 0)
        for row in range(self.row):
            for col in range(self.col):
                data[row][col] = self.data[row][col] == const
        return matrix(data)
#==================================================================================
    def __ne__(self, const):
        if self.shape != other.shape: raise ValueError("Matrix shapes are not matched")     
        data = []
        for i, j in zip(self.data, other.data):
            data.append(i != j)
        return matrix(data)
#==================================================================================
    def __gt__(self, const):
        data = self.__null__(self.shape, key = 0)
        for row in range(self.row):
            for col in range(self.col):
                data[row][col] = self.data[row][col] > const
        return matrix(data)
#==================================================================================
    def __ge__(self, const):
        data = self.__null__(self.shape, key = 0)
        for row in range(self.row):
            for col in range(self.col):
                data[row][col] = self.data[row][col] >= const
        return matrix(data)
#==================================================================================
    def __getitem__(self, index):
        slide = []
        for row in self.data[index[0]]:
            col_ = []
            for col in row:
                col_.append(col[index[1]])
            slide.append(col_)
        return matrix(slide)
#==================================================================================
    def __setitem__(self, index, value):
        self.data[index] = value
#==================================================================================
    def __iter__(self):
        return iter(self.data)

