import numpy as np

def create_symmetric_matrix(n):
    # 先创建一个 n x n 的零矩阵
    mat = np.zeros((n, n), dtype=int)
    
    # 对角线设为1
    np.fill_diagonal(mat, 1)
    
    # 上三角随机 ±1（不包括对角线）
    upper_indices = np.triu_indices(n, k=1)
    mat[upper_indices] = np.random.choice([-1, 1], size=len(upper_indices[0]))
    
    # 下三角赋值为上三角的转置，保证对称
    mat = mat + mat.T - np.diag(np.diag(mat))
    
    return mat

# 举例生成8x8矩阵
# matrix = create_symmetric_matrix(8)
# 定义矩阵
A=np.random.choice([-1, 1], size=8)

M = np.zeros((8, 8))
for i in range (8):
    for j in range (8):
        if i == j:
            M[i][j] = 1
        else:
            M[i][j] = A[i] *A[j]
            M[j][i] = M[i][j]

print(M)
# 计算特征值和特征向量
eigenvalues, eigenvectors = np.linalg.eig(M)

# 找出最大特征值
max_eigenvalue = np.max(eigenvalues)

# print(matrix)
print("最大特征值是:", max_eigenvalue)
